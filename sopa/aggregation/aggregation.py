from __future__ import annotations

import logging

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import TableModel

from .._constants import SopaKeys
from ..io.explorer.utils import str_cell_id
from ..utils import add_spatial_element, get_boundaries, get_spatial_image, to_intrinsic
from . import aggregate_bins, average_channels, count_transcripts

log = logging.getLogger(__name__)


def aggregate(
    sdata: SpatialData,
    average_intensities: bool = True,
    expand_radius_ratio: float = 0,
    min_transcripts: int = 0,
    min_intensity_ratio: float = 0,
    **kwargs: int,
):
    aggr = Aggregator(sdata, **kwargs)

    aggr.compute_table(
        average_intensities=average_intensities,
        expand_radius_ratio=expand_radius_ratio,
        min_transcripts=min_transcripts,
        min_intensity_ratio=min_intensity_ratio,
    )


def overlay_segmentation(
    sdata: SpatialData,
    shapes_key: str,
    gene_column: str | None = None,
    area_ratio_threshold: float = 0.25,
    image_key: str | None = None,
):
    """Overlay a segmentation on top of an existing segmentation

    Args:
        sdata: A `SpatialData` object
        shapes_key: The key of the new shapes to be added
        gene_column: Key of the points dataframe containing the genes names
        area_ratio_threshold: Threshold between 0 and 1. For each original cell overlapping with a new cell, we compute the overlap-area/cell-area, if above the threshold the cell is removed.
        image_key: Optional key of the original image
    """
    average_intensities = False

    if "table" in sdata.tables and SopaKeys.UNS_KEY in sdata.tables["table"].uns:
        sopa_attrs = sdata.tables["table"].uns[SopaKeys.UNS_KEY]

        if sopa_attrs[SopaKeys.UNS_HAS_TRANSCRIPTS]:
            assert gene_column is not None, "Need 'gene_column' argument to count transcripts"
        else:
            gene_column = gene_column
        average_intensities = sopa_attrs[SopaKeys.UNS_HAS_INTENSITIES]

    aggr = Aggregator(sdata, image_key=image_key, shapes_key=shapes_key)
    aggr.overlay_segmentation(
        gene_column=gene_column,
        average_intensities=average_intensities,
        area_ratio_threshold=area_ratio_threshold,
    )


class Aggregator:
    """Perform transcript count and channel averaging over a `SpatialData` object"""

    def __init__(
        self,
        sdata: SpatialData,
        overwrite: bool = True,
        image_key: str | None = None,
        shapes_key: str | None = None,
        bins_key: str | None = None,
    ):
        """
        Args:
            sdata: A `SpatialData` object
            overwrite: If `True`, will overwrite `sdata.table` if already existing
            image_key: Key of `sdata` with the image to be averaged. If only one image, this does not have to be provided
            shapes_key: Key of `sdata` with the shapes corresponding to the cells boundaries
            bins_key: Key of `sdata` with the table corresponding to the bins table of gene counts (e.g., for Visium HD data)
        """
        self.sdata = sdata
        self.overwrite = overwrite

        self.image_key, self.image = get_spatial_image(sdata, image_key, return_key=True)
        self.bins_key = bins_key

        if shapes_key is None:
            self.shapes_key, self.geo_df = get_boundaries(sdata, return_key=True)
        else:
            self.shapes_key = shapes_key
            self.geo_df = self.sdata[shapes_key]

        self.table = None
        self._had_table = False
        if SopaKeys.TABLE in self.sdata.tables:
            table = self.sdata.tables[SopaKeys.TABLE]
            if len(self.geo_df) == table.n_obs:
                log.info("Using existing table for aggregation")
                self.table = table
                self._had_table = True

    def overlay_segmentation(
        self,
        gene_column: str | None = None,
        average_intensities: bool = True,
        area_ratio_threshold: float = 0.25,
    ):
        old_table: AnnData = self.sdata.tables[SopaKeys.TABLE]
        self.sdata.tables[SopaKeys.OLD_TABLE] = old_table
        del self.sdata.tables[SopaKeys.TABLE]

        old_shapes_key = old_table.uns["spatialdata_attrs"]["region"]
        instance_key = old_table.uns["spatialdata_attrs"]["instance_key"]

        if isinstance(old_shapes_key, list):
            assert len(old_shapes_key) == 1, "Can't overlap segmentation on multi-region SpatialData object"
            old_shapes_key = old_shapes_key[0]

        old_geo_df = self.sdata[old_shapes_key]
        geo_df = to_intrinsic(self.sdata, self.geo_df, old_geo_df)

        geo_df.index.name = "index_right"  # to reuse the index name later
        gdf_join = gpd.sjoin(old_geo_df, geo_df)
        gdf_join["geometry_right"] = gdf_join["index_right"].map(lambda i: geo_df.geometry.iloc[i])
        gdf_join["overlap_ratio"] = gdf_join.apply(_overlap_area_ratio, axis=1)
        gdf_join: gpd.GeoDataFrame = gdf_join[gdf_join.overlap_ratio >= area_ratio_threshold]

        table_crop = old_table[~np.isin(old_table.obs[instance_key], gdf_join.index)].copy()
        table_crop.obs[SopaKeys.CELL_OVERLAY_KEY] = False

        self.compute_table(gene_column=gene_column, average_intensities=average_intensities)
        self.table.obs[SopaKeys.CELL_OVERLAY_KEY] = True

        self.table = anndata.concat([table_crop, self.table], uns_merge="first", join="outer", fill_value=0)
        _fillna(self.table.obs)

        self.shapes_key = f"{old_shapes_key}+{self.shapes_key}"
        geo_df_cropped = old_geo_df.loc[~old_geo_df.index.isin(gdf_join.index)]
        self.geo_df = pd.concat([geo_df_cropped, geo_df], join="outer", axis=0)
        self.geo_df.attrs = old_geo_df.attrs

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

        if "spatialdata_attrs" in self.table.uns:
            del self.table.uns["spatialdata_attrs"]

        self.table = TableModel.parse(
            self.table,
            region_key=SopaKeys.REGION_KEY,
            region=self.shapes_key,
            instance_key=SopaKeys.INSTANCE_KEY,
        )

        add_spatial_element(self.sdata, SopaKeys.TABLE, self.table)

    def filter_cells(self, where_filter: np.ndarray):
        log.info(f"Filtering {where_filter.sum()} cells")

        self.geo_df = self.geo_df[~where_filter]

        self.sdata.shapes[self.shapes_key] = self.geo_df

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
                self.table = count_transcripts(self.sdata, gene_column, shapes_key=self.shapes_key)
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

        self.table.uns[SopaKeys.UNS_KEY] = {
            SopaKeys.UNS_HAS_TRANSCRIPTS: does_count,
            SopaKeys.UNS_HAS_INTENSITIES: average_intensities,
        }

        self.add_standardized_table()


def _overlap_area_ratio(row) -> float:
    poly: Polygon = row["geometry"]
    poly_right: Polygon = row["geometry_right"]
    return poly.intersection(poly_right).area / poly.area


def _fillna(df: pd.DataFrame):
    for key in df:
        if df[key].dtype == "category":
            df[key] = df[key].cat.add_categories("NA").fillna("NA")
        else:
            df[key] = df[key].fillna(0)
