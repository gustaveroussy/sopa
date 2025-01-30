import logging

import dask
import geopandas as gpd
import numpy as np
import numpy.ma as ma
import shapely
from dask.diagnostics import ProgressBar
from shapely.geometry import Polygon, box
from spatialdata import SpatialData
from xarray import DataArray

from ..segmentation.shapes import expand_radius, pixel_outer_bounds, rasterize
from ..utils import get_boundaries, get_spatial_image, to_intrinsic

log = logging.getLogger(__name__)

AVAILABLE_MODES = ["average", "min", "max"]


def average_channels(
    sdata: SpatialData, image_key: str = None, shapes_key: str | None = None, expand_radius_ratio: float = 0
) -> np.ndarray:
    log.warning("average_channels is deprecated and will be removed in sopa==2.1.0, use `aggregate_channels` instead")
    return aggregate_channels(sdata, image_key, shapes_key, expand_radius_ratio, mode="average")


def aggregate_channels(
    sdata: SpatialData,
    image_key: str | None = None,
    shapes_key: str | None = None,
    expand_radius_ratio: float = 0,
    mode: str = "average",
) -> np.ndarray:
    """Aggregate the channel intensities per cell (either `"average"`, or take the `"min"` / `"max"`).

    Args:
        sdata: A `SpatialData` object
        image_key: Key of `sdata` containing the image. If only one `images` element, this does not have to be provided.
        shapes_key: Key of `sdata` containing the cell boundaries. If only one `shapes` element, this does not have to be provided.
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate boundary stainings.
        mode: Aggregation mode. One of `"average"`, `"min"`, `"max"`. By default, average the intensity inside the cell mask.

    Returns:
        A numpy `ndarray` of shape `(n_cells, n_channels)`
    """
    assert mode in AVAILABLE_MODES, f"Invalid {mode=}. Available modes are {AVAILABLE_MODES}"

    image = get_spatial_image(sdata, image_key)

    geo_df = get_boundaries(sdata, key=shapes_key)
    geo_df = to_intrinsic(sdata, geo_df, image)
    geo_df = expand_radius(geo_df, expand_radius_ratio)

    return _aggregate_channels_aligned(image, geo_df, mode)


def _aggregate_channels_aligned(image: DataArray, geo_df: gpd.GeoDataFrame | list[Polygon], mode: str) -> np.ndarray:
    """Average channel intensities per cell. The image and cells have to be aligned, i.e. be on the same coordinate system.

    Args:
        image: A `DataArray` of shape `(n_channels, y, x)`
        geo_df: A `GeoDataFrame` whose geometries are cell boundaries (polygons)

    Returns:
        A numpy `ndarray` of shape `(n_cells, n_channels)`
    """
    log.info(f"Aggregating channels intensity over {len(geo_df)} cells with {mode=}")

    cells = geo_df if isinstance(geo_df, list) else list(geo_df.geometry)
    tree = shapely.STRtree(cells)

    n_channels = len(image.coords["c"])
    areas = np.zeros(len(cells))
    if mode == "min":
        aggregation = np.full((len(cells), n_channels), fill_value=np.inf)
    else:
        aggregation = np.zeros((len(cells), n_channels))

    chunk_sizes = image.data.chunks
    offsets_y = np.cumsum(np.pad(chunk_sizes[1], (1, 0), "constant"))
    offsets_x = np.cumsum(np.pad(chunk_sizes[2], (1, 0), "constant"))

    def _average_chunk_inside_cells(chunk, iy, ix):
        ymin, ymax = offsets_y[iy], offsets_y[iy + 1]
        xmin, xmax = offsets_x[ix], offsets_x[ix + 1]

        patch = box(xmin, ymin, xmax, ymax)
        intersections = tree.query(patch, predicate="intersects")

        for index in intersections:
            cell = cells[index]
            bounds = pixel_outer_bounds(cell.bounds)

            sub_image = chunk[
                :,
                max(bounds[1] - ymin, 0) : bounds[3] - ymin,
                max(bounds[0] - xmin, 0) : bounds[2] - xmin,
            ]

            if sub_image.shape[1] == 0 or sub_image.shape[2] == 0:
                continue

            mask = rasterize(cell, sub_image.shape[1:], bounds)

            areas[index] += np.sum(mask)

            if mode == "min":
                masked_image = ma.masked_array(sub_image, 1 - np.repeat(mask[None], n_channels, axis=0))
                aggregation[index] = np.minimum(aggregation[index], masked_image.min(axis=(1, 2)))
            elif mode in ["average", "max"]:
                func = np.sum if mode == "average" else np.max
                values = func(sub_image * mask, axis=(1, 2))

                if mode == "average":
                    aggregation[index] += values
                else:
                    aggregation[index] = np.maximum(aggregation[index], values)

    with ProgressBar():
        tasks = [
            dask.delayed(_average_chunk_inside_cells)(chunk, iy, ix)
            for iy, row in enumerate(image.chunk({"c": -1}).data.to_delayed()[0])
            for ix, chunk in enumerate(row)
        ]
        dask.compute(tasks)

    if mode == "average":
        return aggregation / areas[:, None].clip(1)
    else:
        return aggregation
