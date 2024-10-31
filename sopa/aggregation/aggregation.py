from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix
from spatialdata import SpatialData
from spatialdata.models import PointsModel, TableModel

from .._constants import ATTRS_KEY, SopaAttrs, SopaKeys
from ..io.explorer.utils import str_cell_id
from ..utils import (
    add_spatial_element,
    get_boundaries,
    get_spatial_element,
    get_spatial_image,
)
from . import aggregate_bins, average_channels, count_transcripts

log = logging.getLogger(__name__)


def aggregate(
    sdata: SpatialData,
    aggregate_genes: bool | None = None,
    aggregate_channels: bool = True,
    image_key: str | None = None,
    points_key: str | None = None,
    shapes_key: str | None = None,
    bins_key: str | None = None,
    min_transcripts: int = 0,
    expand_radius_ratio: float | None = None,
    min_intensity_ratio: float = 0.1,
):
    bins_key = bins_key or sdata.attrs.get(SopaAttrs.BINS_TABLE)

    aggr = Aggregator(sdata, image_key=image_key, shapes_key=shapes_key, bins_key=bins_key)

    points_key, gene_column = None, None
    if aggregate_genes or (aggregate_genes is None and bins_key is None and sdata.points):
        points_key, points = get_spatial_element(
            sdata.points, key=points_key or sdata.attrs.get(SopaAttrs.TRANSCRIPTS), return_key=True
        )
        gene_column = points.attrs.get(ATTRS_KEY, {}).get(PointsModel.FEATURE_KEY)
        assert (
            aggregate_genes is None or gene_column is not None
        ), f"No gene column found in sdata['{points_key}'].attrs['{ATTRS_KEY}']['{PointsModel.FEATURE_KEY}']"

    aggr.compute_table(
        gene_column=gene_column,
        average_intensities=aggregate_channels,
        expand_radius_ratio=expand_radius_ratio,
        min_transcripts=min_transcripts,
        min_intensity_ratio=min_intensity_ratio,
        points_key=points_key,
    )


class Aggregator:
    """Perform transcript count and channel averaging over a `SpatialData` object"""

    table: AnnData | None

    def __init__(
        self,
        sdata: SpatialData,
        image_key: str | None = None,
        shapes_key: str | None = None,
        bins_key: str | None = None,
    ):
        """
        Args:
            sdata: A `SpatialData` object
            image_key: Key of `sdata` with the image to be averaged. If only one image, this does not have to be provided
            shapes_key: Key of `sdata` with the shapes corresponding to the cells boundaries
            bins_key: Key of `sdata` with the table corresponding to the bins table of gene counts (e.g., for Visium HD data)
        """
        self.sdata = sdata
        self.bins_key = bins_key

        self.image_key, self.image = get_spatial_image(sdata, image_key, return_key=True)

        if shapes_key is None:
            self.shapes_key, self.geo_df = get_boundaries(sdata, return_key=True)
        else:
            self.shapes_key = shapes_key
            self.geo_df = self.sdata[shapes_key]

        self.table = None
        if SopaKeys.TABLE in self.sdata.tables and self.sdata[SopaKeys.TABLE].n_obs == len(self.geo_df):
            log.info("Using existing table for aggregation")
            self.table = self.sdata[SopaKeys.TABLE]

    def filter_cells(self, where_filter: np.ndarray):
        log.info(f"Filtering {where_filter.sum()} cells")

        self.geo_df = self.geo_df[~where_filter]

        if self.table is not None:
            self.table = self.table[~where_filter]

    def update_table(self, *args, **kwargs):
        log.warning("'update_table' is deprecated, use 'compute_table' instead")
        self.compute_table(*args, **kwargs)

    def compute_table(
        self,
        gene_column: str | None = None,
        average_intensities: bool = True,
        expand_radius_ratio: float | None = None,
        points_key: str | None = None,
        min_transcripts: int = 0,
        min_intensity_ratio: float = 0,
    ):
        """Perform aggregation and update the spatialdata table

        Args:
            gene_column: Column key of the transcript dataframe containing the gene names
            average_intensities: Whether to average the channels intensities inside cells polygons
            expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius` for channels averaging **only**. This help better aggregate boundary stainings
            min_transcripts: Minimum amount of transcript to keep a cell
            min_intensity_ratio: Cells whose mean channel intensity is less than `min_intensity_ratio * quantile_90` will be filtered
        """
        does_count = (
            (self.table is not None and isinstance(self.table.X, csr_matrix))
            or gene_column is not None
            or self.bins_key is not None
        )

        assert (
            average_intensities or does_count
        ), "You must choose at least one aggregation: transcripts or fluorescence intensities"
        assert (
            gene_column is None or self.bins_key is None
        ), "Can't count transcripts and aggregate bins at the same time"

        if gene_column is not None:
            if self.table is not None:
                log.warning("sdata.table is already existing. Transcripts are not count again.")
            else:
                self.table = count_transcripts(
                    self.sdata, gene_column, shapes_key=self.shapes_key, points_key=points_key
                )
        elif self.bins_key is not None:
            self.table = aggregate_bins(
                self.sdata, self.shapes_key, self.bins_key, expand_radius_ratio=expand_radius_ratio or 0.2
            )

        if does_count and min_transcripts > 0:
            self.filter_cells(self.table.X.sum(axis=1) < min_transcripts)

        if average_intensities:
            mean_intensities = average_channels(
                self.sdata,
                image_key=self.image_key,
                shapes_key=self.shapes_key,
                expand_radius_ratio=expand_radius_ratio or 0,
            )

            if min_intensity_ratio > 0:
                means = mean_intensities.mean(axis=1)
                intensity_threshold = min_intensity_ratio * np.quantile(means, 0.9)
                where_filter = means < intensity_threshold
                self.filter_cells(where_filter)
                mean_intensities = mean_intensities[~where_filter]

            if not does_count:
                self.table = AnnData(
                    mean_intensities,
                    dtype=mean_intensities.dtype,
                    var=pd.DataFrame(index=self.image.coords["c"].values.astype(str)),
                    obs=pd.DataFrame(index=self.geo_df.index),
                )
            else:
                self.table.obsm[SopaKeys.INTENSITIES_OBSM] = pd.DataFrame(
                    mean_intensities,
                    columns=self.image.coords["c"].values.astype(str),
                    index=self.table.obs_names,
                )

        self.sdata.shapes[self.shapes_key] = self.geo_df

        self.table.uns[SopaKeys.UNS_KEY] = {
            SopaKeys.UNS_HAS_TRANSCRIPTS: does_count,
            SopaKeys.UNS_HAS_INTENSITIES: average_intensities,
        }

        self.add_standardized_table()

    def add_standardized_table(self):
        self.table.obs_names = list(map(str_cell_id, range(self.table.n_obs)))
        self.geo_df.index = list(self.table.obs_names)

        add_spatial_element(self.sdata, self.shapes_key, self.geo_df)

        self.table.obsm["spatial"] = np.array([[centroid.x, centroid.y] for centroid in self.geo_df.centroid])
        self.table.obs[SopaKeys.REGION_KEY] = pd.Series(self.shapes_key, index=self.table.obs_names, dtype="category")
        self.table.obs[SopaKeys.SLIDE_KEY] = pd.Series(self.image_key, index=self.table.obs_names, dtype="category")
        self.table.obs[SopaKeys.INSTANCE_KEY] = self.geo_df.index
        self.table.obs[SopaKeys.AREA_OBS] = self.geo_df.area.values

        if ATTRS_KEY in self.table.uns:
            del self.table.uns[ATTRS_KEY]

        self.table = TableModel.parse(
            self.table,
            region_key=SopaKeys.REGION_KEY,
            region=self.shapes_key,
            instance_key=SopaKeys.INSTANCE_KEY,
        )

        add_spatial_element(self.sdata, SopaKeys.TABLE, self.table)
