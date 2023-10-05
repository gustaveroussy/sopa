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
from shapely.geometry import Point, Polygon, box
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import TableModel

from .._constants import SopaKeys
from .._sdata import (
    get_boundaries,
    get_intrinsic_cs,
    get_item,
    get_key,
    get_spatial_image,
    to_intrinsic,
)
from . import shapes

log = logging.getLogger(__name__)


def _average_channels(image: SpatialImage, cells: list[Polygon]):
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
    sdata: SpatialData, image: SpatialImage, geo_df: gpd.GeoDataFrame
) -> np.ndarray:
    geo_df = to_intrinsic(sdata, geo_df, image)
    cells = geo_df.geometry

    log.info(f"Averaging intensities over {len(cells)} cells")
    return _average_channels(image, cells)


def _get_cell_id(polygons: list[Polygon], partition: pd.DataFrame) -> np.ndarray:
    points = partition[["x", "y"]].apply(Point, axis=1)
    tree = shapely.STRtree(points)

    polygon_indices, points_indices = tree.query(polygons, predicate="contains")

    unique_values, indices = np.unique(points_indices, return_index=True)
    cell_id = np.zeros(len(points), dtype=int)
    cell_id[unique_values] = 1 + polygon_indices[indices]

    return cell_id


def map_transcript_to_cell(
    sdata: SpatialData,
    cell_key: str,
    df: dd.DataFrame | None = None,
    polygons: list[Polygon] | None = None,
):
    if df is None:
        points_key, df = get_item(sdata, "points")

    if polygons is None:
        geo_df = get_boundaries(sdata)

    polygons = to_intrinsic(sdata, geo_df, points_key).geometry

    get_cell_id = partial(_get_cell_id, polygons)

    if isinstance(df, dd.DataFrame):
        df[cell_key] = df.map_partitions(get_cell_id)
    else:
        raise ValueError(f"Invalid dataframe type: {type(df)}")


def aggregate(
    sdata: SpatialData, gene_column: str | None, intensity_mean: bool = True, overwrite: bool = True
):
    image_key, image = get_spatial_image(sdata)
    shapes_key, geo_df = get_boundaries(sdata, return_key=True)

    table = table if sdata.table is None else sdata.table

    assert (
        intensity_mean or gene_column is not None or table is not None
    ), f"You must choose at least one aggregation: transcripts or fluorescence intensities"

    if gene_column is not None:
        points_key = get_key(sdata, "points")

        table = sdata.aggregate(
            values=points_key,
            by=shapes_key,
            value_key=gene_column,
            agg_func="count",
            target_coordinate_system=get_intrinsic_cs(sdata, points_key),
        ).table

    if intensity_mean:
        mean_intensities = average_channels(sdata, image, geo_df)

    if table is None:
        table = AnnData(
            mean_intensities,
            dtype=mean_intensities.dtype,
            var=pd.DataFrame(index=image.c),
            obs=pd.DataFrame(index=geo_df.index),
        )
    elif intensity_mean:
        table.obsm[SopaKeys.INTENSITIES_OBSM] = pd.DataFrame(
            mean_intensities, columns=image.coords["c"].values, index=table.obs_names
        )

    table.obsm["spatial"] = np.array([[centroid.x, centroid.y] for centroid in geo_df.centroid])
    table.obs[SopaKeys.REGION_KEY] = pd.Series(shapes_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.SLIDE_KEY] = pd.Series(image_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.INSTANCE_KEY] = geo_df.index

    if "spatialdata_attrs" in table.uns:
        del table.uns["spatialdata_attrs"]

    table = TableModel.parse(
        table,
        region_key=SopaKeys.REGION_KEY,
        region=shapes_key,
        instance_key=SopaKeys.INSTANCE_KEY,
    )

    if sdata.table is not None and overwrite:
        del sdata.table

    sdata.table = table
