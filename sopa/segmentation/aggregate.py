import logging
from functools import partial

import dask.array as da
import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from anndata import AnnData
from dask.diagnostics import ProgressBar
from scipy.sparse import coo_matrix, csr_matrix
from shapely.geometry import Polygon, box
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import TableModel

import sopa

from .._constants import SopaKeys
from .._sdata import (
    get_boundaries,
    get_element,
    get_item,
    get_spatial_image,
    to_intrinsic,
)
from ..io.explorer.utils import str_cell_id
from . import shapes

log = logging.getLogger(__name__)


class Aggregator:
    """Perform transcript count and channel averaging over a `SpatialData` object"""

    def __init__(
        self,
        sdata: SpatialData,
        overwrite: bool = True,
        image_key: str | None = None,
        shapes_key: str | None = None,
    ):
        """
        Args:
            sdata: A `SpatialData` object
            overwrite: If `True`, will overwrite `sdata.table` if already existing
            image_key: Key of `sdata` with the image to be averaged. If only one image, this does not have to be provided.
            shapes_key: Key of `sdata` with the shapes corresponding to the cells boundaries
        """
        self.sdata = sdata
        self.overwrite = overwrite

        self.image_key, self.image = get_spatial_image(sdata, image_key, return_key=True)

        if shapes_key is None:
            self.shapes_key, self.geo_df = get_boundaries(sdata, return_key=True)
        else:
            self.shapes_key = shapes_key
            self.geo_df = self.sdata[shapes_key]

        if sdata.table is not None and len(self.geo_df) != sdata.table.n_obs:
            log.warn(
                f"Table already existing with {sdata.table.n_obs} obs, but aggregating on {len(self.geo_df)} cells. Deleting table."
            )
            del sdata.table

        self.table = sdata.table

    def standardize_table(self):
        self.table.obs_names = list(map(str_cell_id, range(self.table.n_obs)))

        self.table.obsm["spatial"] = np.array(
            [[centroid.x, centroid.y] for centroid in self.geo_df.centroid]
        )
        self.table.obs[SopaKeys.REGION_KEY] = pd.Series(
            self.shapes_key, index=self.table.obs_names, dtype="category"
        )
        self.table.obs[SopaKeys.SLIDE_KEY] = pd.Series(
            self.image_key, index=self.table.obs_names, dtype="category"
        )
        self.table.obs[SopaKeys.INSTANCE_KEY] = self.geo_df.index

        if "spatialdata_attrs" in self.table.uns:
            del self.table.uns["spatialdata_attrs"]

        self.table = TableModel.parse(
            self.table,
            region_key=SopaKeys.REGION_KEY,
            region=self.shapes_key,
            instance_key=SopaKeys.INSTANCE_KEY,
        )

    def filter_cells(self, where_filter: np.ndarray):
        log.info(f"Filtering {where_filter.sum()} cells")

        self.geo_df = self.geo_df[~where_filter]
        self.sdata.add_shapes(self.shapes_key, self.geo_df, overwrite=True)

        if self.table is not None:
            self.table = self.table[~where_filter]

    def save_table(self):
        if self.sdata.table is not None and self.overwrite:
            del self.sdata.table
        self.sdata.table = self.table

    def update_table(
        self,
        gene_column: str | None,
        average_intensities: bool,
        expand_radius_ratio: float,
        min_transcripts: int,
        min_intensity_ratio: float,
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
            self.table is not None and isinstance(self.table.X, csr_matrix)
        ) or gene_column is not None

        assert (
            average_intensities or does_count
        ), f"You must choose at least one aggregation: transcripts or fluorescence intensities"

        if gene_column is not None:
            if self.table is not None:
                log.warn("sdata.table is already existing. Transcripts are not count again.")
            else:
                self.table = count_transcripts(self.sdata, gene_column, shapes_key=self.shapes_key)

        if does_count and min_transcripts > 0:
            self.filter_cells(self.table.X.sum(axis=1) < min_transcripts)

        if average_intensities:
            mean_intensities = average_channels(
                self.sdata,
                image_key=self.image_key,
                shapes_key=self.shapes_key,
                expand_radius_ratio=expand_radius_ratio,
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
            "version": sopa.__version__,
            SopaKeys.UNS_HAS_TRANSCRIPTS: does_count,
            SopaKeys.UNS_HAS_INTENSITIES: average_intensities,
        }

        self.standardize_table()
        self.save_table()


def average_channels(
    sdata: SpatialData,
    image_key: str = None,
    shapes_key: str = None,
    expand_radius_ratio: float = 0,
) -> np.ndarray:
    """Average channel intensities per cell.

    Args:
        sdata: A `SpatialData` object
        image_key: Key of `sdata` containing the image. If only one `images` element, this does not have to be provided.
        shapes_key: Key of `sdata` containing the cell boundaries. If only one `shapes` element, this does not have to be provided.
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate boundary stainings.

    Returns:
        A numpy `ndarray` of shape `(n_cells, n_channels)`
    """
    image = get_spatial_image(sdata, image_key)

    geo_df = get_element(sdata, "shapes", shapes_key)
    geo_df = to_intrinsic(sdata, geo_df, image)

    expand_radius = expand_radius_ratio * np.mean(np.sqrt(geo_df.area / np.pi))

    if expand_radius > 0:
        geo_df = geo_df.buffer(expand_radius)

    log.info(
        f"Averaging channels intensity over {len(geo_df)} cells with expansion {expand_radius}"
    )
    return _average_channels_aligned(image, geo_df)


def _average_channels_aligned(
    image: SpatialImage, geo_df: gpd.GeoDataFrame | list[Polygon]
) -> np.ndarray:
    """Average channel intensities per cell. The image and cells have to be aligned, i.e. be on the same coordinate system.

    Args:
        image: A `SpatialImage` of shape `(n_channels, y, x)`
        geo_df: A `GeoDataFrame` whose geometries are cell boundaries (polygons)

    Returns:
        A numpy `ndarray` of shape `(n_cells, n_channels)`
    """
    cells = geo_df if isinstance(geo_df, list) else list(geo_df.geometry)
    tree = shapely.STRtree(cells)

    intensities = np.zeros((len(cells), len(image.coords["c"])))
    areas = np.zeros(len(cells))

    def func(chunk, block_info=None):
        if block_info is not None:
            (ymin, ymax), (xmin, xmax) = block_info[0]["array-location"][1:]
            patch = box(xmin, ymin, xmax, ymax)
            intersections = tree.query(patch, predicate="intersects")

            for index in intersections:
                cell = cells[index]
                bounds = shapes.pixel_outer_bounds(cell.bounds)

                sub_image = chunk[
                    :,
                    max(bounds[1] - ymin, 0) : bounds[3] - ymin,
                    max(bounds[0] - xmin, 0) : bounds[2] - xmin,
                ]

                if sub_image.shape[1] == 0 or sub_image.shape[2] == 0:
                    continue

                mask = shapes.rasterize(cell, sub_image.shape[1:], bounds)

                intensities[index] += np.sum(sub_image * mask, axis=(1, 2))
                areas[index] += np.sum(mask)
        return da.zeros(chunk.shape[1:])

    with ProgressBar():
        image.data.map_blocks(func, drop_axis=0).compute()

    return intensities / areas[:, None].clip(1)


def count_transcripts(
    sdata: SpatialData,
    gene_column: str,
    shapes_key: str = None,
    points_key: str = None,
    geo_df: gpd.GeoDataFrame = None,
) -> AnnData:
    """Counts transcripts per cell.

    Args:
        sdata: A `SpatialData` object
        gene_column: Column of the transcript dataframe containing the gene names
        shapes_key: Key of `sdata` containing the cell boundaries. If only one `shapes` element, this does not have to be provided.
        points_key: Key of `sdata` containing the transcripts. If only one `points` element, this does not have to be provided.
        geo_df: If the cell boundaries are not yet in `sdata`, a `GeoDataFrame` can be directly provided for cell boundaries

    Returns:
        An `AnnData` object of shape `(n_cells, n_genes)` with the counts per cell
    """
    points_key, points = get_item(sdata, "points", points_key)

    if geo_df is None:
        geo_df = get_element(sdata, "shapes", shapes_key)
        geo_df = to_intrinsic(sdata, geo_df, points_key)

    log.info(f"Aggregating transcripts over {len(geo_df)} cells")
    return _count_transcripts_aligned(geo_df, points, gene_column)


def _count_transcripts_aligned(
    geo_df: gpd.GeoDataFrame, points: dd.DataFrame, value_key: str
) -> AnnData:
    """Count transcripts per cell. The cells and points have to be aligned (i.e., in the same coordinate system)

    Args:
        geo_df: Cells geometries
        points: Transcripts dataframe
        value_key: Key of `points` containing the genes names

    Returns:
        An `AnnData` object of shape `(n_cells, n_genes)` with the counts per cell
    """
    points[value_key] = points[value_key].astype("category").cat.as_known()
    gene_names = points[value_key].cat.categories.astype(str)

    X = coo_matrix((len(geo_df), len(gene_names)), dtype=int)
    adata = AnnData(X=X, var=pd.DataFrame(index=gene_names))
    adata.obs_names = geo_df.index

    geo_df = geo_df.reset_index()

    with ProgressBar():
        points.map_partitions(
            partial(_add_coo, adata, geo_df, gene_column=value_key, gene_names=gene_names),
            meta=(),
        ).compute()

    adata.X = adata.X.tocsr()
    return adata


def _add_coo(
    adata: AnnData,
    geo_df: gpd.GeoDataFrame,
    partition: pd.DataFrame,
    gene_column: list[str],
    gene_names: list[str],
) -> None:
    points_gdf = gpd.GeoDataFrame(
        partition, geometry=gpd.points_from_xy(partition["x"], partition["y"])
    )
    joined = geo_df.sjoin(points_gdf)
    cells_indices, column_indices = joined.index, joined[gene_column].cat.codes

    adata.X += coo_matrix(
        (np.full(len(cells_indices), 1), (cells_indices, column_indices)),
        shape=(len(geo_df), len(gene_names)),
    )
