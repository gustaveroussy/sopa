import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix
from spatialdata import SpatialData
from spatialdata.models import TableModel

from .._constants import ATTRS_KEY, SopaAttrs, SopaKeys
from ..io.explorer.utils import str_cell_id
from ..utils import (
    add_spatial_element,
    get_boundaries,
    get_spatial_element,
    get_spatial_image,
)
from . import aggregate_bins
from . import aggregate_channels as _aggregate_channels
from . import count_transcripts

log = logging.getLogger(__name__)


def aggregate(
    sdata: SpatialData,
    aggregate_genes: bool | None = None,
    aggregate_channels: bool = True,
    image_key: str | None = None,
    points_key: str | None = None,
    gene_column: str | None = None,
    shapes_key: str | None = None,
    bins_key: str | None = None,
    expand_radius_ratio: float | None = None,
    min_transcripts: int = 0,
    min_intensity_ratio: float = 0.1,
    key_added: str | None = "table",
):
    """Aggregate gene counts and/or channel intensities over a `SpatialData` object to create an `AnnData` table (saved in `sdata["table"]`).

    !!! info
        The main arguments are `sdata`, `aggregate_genes`, and `aggregate_channels`. The rest of the arguments are optional and will be inferred from the data if not provided.

        - If channels are aggregated and not genes, then `sdata['table'].X` will contain the mean channel intensities per cell.
        - If genes are aggregated and not channels, then `sdata['table'].X` will contain the gene counts per cell.
        - If both genes and channels are aggregated, then `sdata['table'].X` will contain the gene counts per cell and `sdata['table'].obsm['intensities']` will contain the mean channel intensities per cell.

    Args:
        sdata: A `SpatialData` object
        aggregate_genes: Whether to aggregate gene counts. If None, it will be inferred from the data.
        aggregate_channels: Whether to aggregate channel intensities inside cells.
        image_key: Key of `sdata` with the image channels to be averaged. By default, uses the segmentation image.
        points_key: Key of `sdata` with the points dataframe representing the transcripts.
        gene_column: Key of `sdata[points_key]` with the gene names.
        shapes_key: Key of `sdata` with the shapes corresponding to the cells boundaries.
        bins_key: Key of `sdata` with the table corresponding to the bin-by-gene table of gene counts (e.g., for Visium HD data).
        expand_radius_ratio: Ratio to expand the cells polygons for channels averaging. For instance, a ratio of 0.5 expands the shape radius by 50%. If `None` (default), use 1 if we aggregate bins data, and 0 otherwise.
        min_transcripts: Min number of transcripts to keep a cell.
        min_intensity_ratio: Min ratio of the 90th quantile of the mean channel intensity to keep a cell.
        key_added: Key to save the table in `sdata.tables`. If `None`, it will be `f"{shapes_key}_table"`.
    """
    assert points_key is None or bins_key is None, "Provide either `points_key` or `bins_key`, not both."

    if points_key is None:
        bins_key = bins_key or sdata.attrs.get(SopaAttrs.BINS_TABLE)

    if (bins_key is None) and (aggregate_genes or (aggregate_genes is None and sdata.points)):
        assert (
            sdata.points
        ), "No points in the SpatialData object. You must have points, or set the `bins_key` argument (for VisiumHD-like data)."

        points_key, _ = get_spatial_element(
            sdata.points, key=points_key or sdata.attrs.get(SopaAttrs.TRANSCRIPTS), return_key=True
        )

    aggr = Aggregator(sdata, image_key=image_key, shapes_key=shapes_key, bins_key=bins_key, points_key=points_key)

    if key_added is None:
        key_added = f"{aggr.shapes_key}_{SopaKeys.TABLE}"
        log.info(f"key_added is None, saving the table as '{key_added}' by default.")

    aggr.compute_table(
        aggregate_genes=aggregate_genes,
        aggregate_channels=aggregate_channels,
        expand_radius_ratio=expand_radius_ratio,
        min_transcripts=min_transcripts,
        min_intensity_ratio=min_intensity_ratio,
        gene_column=gene_column,
        key_added=key_added,
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
        points_key: str | None = None,
    ):
        """
        Args:
            sdata: A `SpatialData` object
            image_key: Key of `sdata` with the image to be averaged. If only one image, this does not have to be provided
            shapes_key: Key of `sdata` with the shapes corresponding to the cells boundaries
            bins_key: Key of `sdata` with the table corresponding to the bins table of gene counts (e.g., for Visium HD data)
            points_key: Key of `sdata` with the points dataframe representing the transcripts
        """
        self.sdata = sdata
        self.bins_key = bins_key
        self.points_key = points_key
        self.shapes_key, self.geo_df = get_boundaries(sdata, return_key=True, key=shapes_key)
        self.table = None

        if not sdata.images:
            self.image_key, self.image = "None", None
        else:
            self.image_key, self.image = get_spatial_image(sdata, image_key, return_key=True)

    def already_has_valid_table(self, key_added: str) -> bool:
        if (key_added not in self.sdata.tables) or (self.sdata[key_added].n_obs != len(self.geo_df)):
            return False

        log.info("Found existing table, transcripts are not count again.")
        self.table = self.sdata[key_added]

        return True

    def filter_cells(self, where_filter: np.ndarray):
        log.info(f"Filtering {where_filter.sum()} cells")

        self.geo_df = self.geo_df[~where_filter]
        self.sdata.shapes[self.shapes_key] = self.geo_df

        if self.table is not None:
            self.table = self.table[~where_filter]

    def update_table(self, *args, **kwargs):
        log.warning("'update_table' is deprecated and will be removed in sopa==2.1.0, use 'compute_table' instead")
        self.compute_table(*args, **kwargs)

    def compute_table(
        self,
        aggregate_genes: bool | None = None,
        aggregate_channels: bool = True,
        gene_column: str | None = None,
        expand_radius_ratio: float | None = None,
        min_transcripts: int = 0,
        min_intensity_ratio: float = 0,
        average_intensities: bool | None = None,  # deprecated argument
        points_key: str | None = None,  # deprecated argument
        key_added: str = SopaKeys.TABLE,
    ):
        aggregate_genes, aggregate_channels = self._legacy_arguments(
            points_key, gene_column, aggregate_genes, aggregate_channels, average_intensities
        )

        if aggregate_channels and self.image is None:
            log.warning("No image found to aggregate channels. Use `aggregate_channels=False`")
            aggregate_channels = False

        assert (
            aggregate_genes or aggregate_channels
        ), "At least one of `aggregate_genes` or `aggregate_channels` must be True"

        if aggregate_genes:
            if self.bins_key is not None:
                self.table = aggregate_bins(
                    self.sdata, self.shapes_key, self.bins_key, expand_radius_ratio=expand_radius_ratio or 1
                )
            elif not self.already_has_valid_table(key_added):
                self.table = count_transcripts(
                    self.sdata, gene_column, shapes_key=self.shapes_key, points_key=points_key
                )

            if min_transcripts > 0:
                self.filter_cells(self.table.X.sum(axis=1) < min_transcripts)

        if aggregate_channels:
            mean_intensities = _aggregate_channels(
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

            if aggregate_genes:
                self.table.obsm[SopaKeys.INTENSITIES_OBSM] = pd.DataFrame(
                    mean_intensities,
                    columns=self.image.coords["c"].values.astype(str),
                    index=self.table.obs_names,
                )
            else:
                self.table = AnnData(
                    mean_intensities,
                    var=pd.DataFrame(index=self.image.coords["c"].values.astype(str)),
                    obs=pd.DataFrame(index=self.geo_df.index),
                )

        self.sdata.shapes[self.shapes_key] = self.geo_df

        self.table.uns[SopaKeys.UNS_KEY] = {
            SopaKeys.UNS_HAS_TRANSCRIPTS: aggregate_genes,
            SopaKeys.UNS_HAS_INTENSITIES: aggregate_channels,
        }

        self.add_standardized_table(key_added)

    def _legacy_arguments(
        self,
        points_key: str | None,
        gene_column: str | None,
        aggregate_genes: bool | None,
        aggregate_channels: bool,
        average_intensities: bool | None,
    ) -> tuple[bool, bool]:
        if points_key is not None:
            log.warning(
                "`points_key` in `compute_table` is deprecated and will be removed in sopa==2.1.0, provide it in the constructor instead"
            )
            self.points_key = points_key

        if aggregate_genes is None:
            aggregate_genes = (
                (self.table is not None and isinstance(self.table.X, csr_matrix))
                or gene_column is not None
                or self.bins_key is not None
                or self.points_key is not None
            )

        if average_intensities is not None:
            log.warning(
                "`average_intensities` is deprecated and will be removed in sopa==2.1.0, use `aggregate_channels` instead"
            )
            return aggregate_genes, average_intensities

        return aggregate_genes, aggregate_channels

    def add_standardized_table(self, key_added: str):
        add_standardized_table(
            self.sdata, self.table, self.geo_df, self.shapes_key, key_added, image_key=self.image_key
        )


def add_standardized_table(
    sdata: SpatialData,
    table: AnnData,
    geo_df: gpd.GeoDataFrame,
    shapes_key: str,
    table_key: str,
    image_key: str | None = None,
):
    table.obs_names = list(map(str_cell_id, range(table.n_obs)))
    geo_df.index = list(table.obs_names)

    add_spatial_element(sdata, shapes_key, geo_df)

    table.obsm["spatial"] = np.array([[centroid.x, centroid.y] for centroid in geo_df.centroid])
    table.obs[SopaKeys.REGION_KEY] = pd.Series(shapes_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.SLIDE_KEY] = pd.Series(image_key or "None", index=table.obs_names, dtype="category")
    table.obs[SopaKeys.INSTANCE_KEY] = geo_df.index
    table.obs[SopaKeys.AREA_OBS] = geo_df.area.values

    if ATTRS_KEY in table.uns:
        del table.uns[ATTRS_KEY]

    table = TableModel.parse(
        table,
        region_key=SopaKeys.REGION_KEY,
        region=shapes_key,
        instance_key=SopaKeys.INSTANCE_KEY,
    )

    add_spatial_element(sdata, table_key, table)
