import logging
from functools import partial

import dask.array as da
import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import shapely.affinity
from anndata import AnnData
from dask.diagnostics import ProgressBar
from scipy.sparse import coo_matrix
from shapely.geometry import Point, Polygon, box
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
from . import shapes

log = logging.getLogger(__name__)


class Aggregator:
    def __init__(self, sdata: SpatialData, overwrite: bool = True):
        self.sdata = sdata
        self.overwrite = overwrite

        self.image_key, self.image = get_spatial_image(sdata, return_key=True)
        self.shapes_key, self.geo_df = get_boundaries(sdata, return_key=True)

        self.table = sdata.table

    def standardize_table(self):
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

    def update_table(
        self,
        gene_column: str | None,
        average_intensities: bool,
        min_transcripts: int,
        min_intensity_ratio: float,
    ):
        does_count = self.table is not None or gene_column is not None

        assert (
            average_intensities or does_count
        ), f"You must choose at least one aggregation: transcripts or fluorescence intensities"

        if gene_column is not None:
            if self.table is not None:
                log.warn("sdata.table is already existing. Transcripts are not count.")
            else:
                self.table = count_transcripts(self.sdata, gene_column, shapes_key=self.shapes_key)

        if does_count and min_transcripts > 0:
            self.filter_cells(self.table.X.sum(axis=1) < min_transcripts)

        if average_intensities:
            mean_intensities = average_channels(self.sdata, shapes_key=self.shapes_key)

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
                    var=pd.DataFrame(index=self.image.c),
                    obs=pd.DataFrame(index=self.geo_df.index),
                )
            else:
                self.table.obsm[SopaKeys.INTENSITIES_OBSM] = pd.DataFrame(
                    mean_intensities,
                    columns=self.image.coords["c"].values,
                    index=self.table.obs_names,
                )

        self.table.uns["sopa_attrs"] = {
            "version": sopa.__version__,
            "transcripts": does_count,
            "intensities": average_intensities,
        }

        self.standardize_table()

        if self.sdata.table is not None and self.overwrite:
            del self.sdata.table

        self.sdata.table = self.table


def _average_channels_geometries(image: SpatialImage, cells: list[Polygon]):
    tree = shapely.STRtree(cells)

    intensities = np.zeros((len(cells), len(image.coords["c"])))
    areas = np.zeros(len(cells))

    def func(x, block_info=None):
        if block_info is not None:
            (ymin, ymax), (xmin, xmax) = block_info[0]["array-location"][1:]
            patch = box(xmin, ymin, xmax, ymax)
            intersections = tree.query(patch, predicate="intersects")

            for index in intersections:
                cell = cells[index]
                bounds = shapes.pixel_outer_bounds(cell.bounds)

                sub_image = x[
                    :,
                    max(bounds[1] - ymin, 0) : bounds[3] - ymin,
                    max(bounds[0] - xmin, 0) : bounds[2] - xmin,
                ]

                if sub_image.shape[1] == 0 or sub_image.shape[2] == 0:
                    continue

                mask = shapes.rasterize(cell, sub_image.shape[1:], bounds)

                intensities[index] += np.sum(sub_image * mask, axis=(1, 2))
                areas[index] += np.sum(mask)
        return da.zeros_like(x)

    image.data.rechunk({0: -1}).map_blocks(func).compute()

    return intensities / areas[:, None].clip(1)


def average_channels(
    sdata: SpatialData, image_key: str = None, shapes_key: str = None
) -> np.ndarray:
    image = get_spatial_image(sdata, image_key)

    geo_df = get_element(sdata, "shapes", shapes_key)
    geo_df = to_intrinsic(sdata, geo_df, image)

    cells = list(geo_df.geometry)

    log.info(f"Averaging channels intensity over {len(cells)} cells")
    return _average_channels_geometries(image, cells)


def _count_transcripts_geometries(
    polygons: list[Polygon], points: dd.DataFrame, value_key: str
) -> AnnData:
    points[value_key] = points[value_key].astype("category").cat.as_known()
    gene_names = points[value_key].cat.categories

    X = coo_matrix((len(polygons), len(gene_names)), dtype=int)
    adata = AnnData(X=X, var=pd.DataFrame(index=gene_names))

    with ProgressBar():
        points.map_partitions(
            partial(_add_coo, adata, polygons, gene_column=value_key, gene_names=gene_names),
            meta=(),
        ).compute()

    adata.X = adata.X.tocsr()
    return adata


def count_transcripts(
    sdata: SpatialData,
    gene_column: str,
    shapes_key: str = None,
    points_key: str = None,
    geo_df: gpd.GeoDataFrame = None,
) -> AnnData:
    points_key, points = get_item(sdata, "points", points_key)

    if geo_df is None:
        geo_df = get_element(sdata, "shapes", shapes_key)
        geo_df = to_intrinsic(sdata, geo_df, points_key)

    cells = list(geo_df.geometry)

    log.info(f"Aggregating transcripts over {len(cells)} cells")
    return _count_transcripts_geometries(cells, points, gene_column)


def _add_coo(
    adata: AnnData,
    cells: list[Polygon],
    partition: pd.DataFrame,
    gene_column: list[str],
    gene_names: list[str],
) -> None:
    cells_indices, points_indices = _get_cell_id(cells, partition, return_query=True)

    column_indices = partition.iloc[points_indices][gene_column].cat.codes.values

    adata.X += coo_matrix(
        (np.full(len(cells_indices), 1), (cells_indices, column_indices)),
        shape=(len(cells), len(gene_names)),
    )


def _get_cell_id(
    polygons: list[Polygon], partition: pd.DataFrame, return_query: bool = False, na_cells: int = 0
) -> np.ndarray:
    points = partition[["x", "y"]].apply(Point, axis=1)
    tree = shapely.STRtree(points)

    polygon_indices, points_indices = tree.query(polygons, predicate="intersects")
    if return_query:
        return polygon_indices, points_indices

    unique_values, indices = np.unique(points_indices, return_index=True)
    cell_id = np.full(len(points), na_cells, dtype=int)

    cell_id[unique_values] = polygon_indices[indices] + 1 + na_cells

    return cell_id


def _map_transcript_to_cell(
    sdata: SpatialData,
    cell_key: str,
    df: dd.DataFrame | None = None,
    geo_df: gpd.GeoDataFrame | None = None,
):
    if df is None:
        df = get_item(sdata, "points")

    if geo_df is None:
        geo_df = get_boundaries(sdata)

    polygons = to_intrinsic(sdata, geo_df, df).geometry

    get_cell_id = partial(_get_cell_id, polygons)

    if isinstance(df, dd.DataFrame):
        df[cell_key] = df.map_partitions(get_cell_id)
    else:
        raise ValueError(f"Invalid dataframe type: {type(df)}")
