import logging
from math import ceil, floor

import numpy as np
import shapely
import shapely.affinity
import xarray as xr
from shapely.geometry import MultiPolygon, Polygon, box

log = logging.getLogger(__name__)


def solve_conflicts(
    polygons: list[Polygon],
    threshold: float = 0.5,
    patch_indices: np.ndarray | None = None,
    return_indices: bool = False,
) -> np.ndarray[Polygon]:
    n_polygons = len(polygons)
    resolved_indices = np.arange(n_polygons)

    tree = shapely.STRtree(polygons)
    conflicts = tree.query(polygons, predicate="intersects")

    if patch_indices is not None:
        conflicts = conflicts[:, patch_indices[conflicts[0]] != patch_indices[conflicts[1]]].T
    else:
        conflicts = conflicts[:, conflicts[0] != conflicts[1]].T

    for i1, i2 in conflicts:
        resolved_i1, resolved_i2 = resolved_indices[i1], resolved_indices[i2]
        poly1, poly2 = polygons[resolved_i1], polygons[resolved_i2]

        intersection = poly1.intersection(poly2).area
        if intersection / min(poly1.area, poly2.area) >= threshold:
            resolved_indices[np.isin(resolved_indices, [resolved_i1, resolved_i2])] = len(polygons)
            polygons.append(poly1.union(poly2))

    unique_indices = np.unique(resolved_indices)
    unique_polygons = np.array(polygons)[unique_indices]

    if return_indices:
        return unique_polygons, np.where(unique_indices < n_polygons, unique_indices, -1)

    return unique_polygons


def expand(polygons: list[Polygon], expand_radius: float) -> list[Polygon]:
    return [polygon.buffer(expand_radius) for polygon in polygons]


def smooth(poly: Polygon, dist: float, tolerance: float) -> Polygon:
    return poly.buffer(-dist).buffer(-2 * dist).buffer(-dist).simplify(tolerance)


def _find_contours(cell_mask: np.ndarray) -> list[Polygon]:
    import cv2

    contours, _ = cv2.findContours(cell_mask, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    return [Polygon(c[:, 0, :]) for c in contours if c.shape[0] >= 4]


def _geometrize_cell(
    mask: np.ndarray, cell_id: int, smooth_radius: int, tolerance: float
) -> Polygon | None:
    """Transform one cell id to a polygon based on a mask array. Inspired from https://github.com/Vizgen/vizgen-postprocessing

    Args:
        mask: Numpy array with values == cell_id where the cell is
        cell_id: ID of the cell to geometrize
        smooth_radius: radius used to smooth the cell polygon
        tolerance: tolerance used to simplify the cell polygon

    Returns:
        Shapely polygon representing the cell, or `None` if the cell was empty after smoothing
    """
    polygons = _find_contours((mask == cell_id).astype("uint8"))
    polygons = [polygon for polygon in polygons if not polygon.buffer(-smooth_radius).is_empty]
    polygon = MultiPolygon(polygons)
    polygon = polygon.buffer(smooth_radius).buffer(-smooth_radius).simplify(tolerance)

    if polygon.is_empty:
        return None

    if isinstance(polygon, Polygon):
        return polygon

    log.warn(
        f"""Geometry index {cell_id} is composed of {len(polygon.geoms)} polygons of areas: {[p.area for p in polygon.geoms]}. Only the polygon corresponding to the largest area will be kept"""
    )

    return max(polygon.geoms, key=lambda polygon: polygon.area)


def geometrize(mask: np.ndarray, smooth_radius: int = 5, tolerance: float = 2) -> list[Polygon]:
    polygons = [
        _geometrize_cell(mask, cell_id, smooth_radius, tolerance)
        for cell_id in range(1, mask.max() + 1)
    ]
    return [polygon for polygon in polygons if polygon is not None]


def pixel_outer_bounds(bounds: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return [floor(bounds[0]), floor(bounds[1]), ceil(bounds[2]) + 1, ceil(bounds[3]) + 1]


def rasterize(
    poly: Polygon, shape: tuple[int, int], xy_min: tuple[int, int] = [0, 0]
) -> np.ndarray:
    """Transform a polygon into a numpy array with value 1 where the polygon touches a pixel, else 0.

    Args:
        poly: Polygon to rasterize.
        shape: Image shape as a tuple (y, x).
        xy_min: Tuple containing the origin of the image [x0, y0].

    Returns:
        The mask array.
    """
    import cv2

    xmin, ymin, xmax, ymax = [xy_min[0], xy_min[1], xy_min[0] + shape[1], xy_min[1] + shape[0]]

    new_poly = shapely.affinity.translate(poly, -xmin, -ymin)
    coords = np.array(new_poly.boundary.coords)[None, :].astype(np.int32)
    return cv2.fillPoly(np.zeros((ymax - ymin, xmax - xmin), dtype=np.int8), coords, color=1)


def average(xarr: xr.DataArray, cells: list[Polygon]) -> np.ndarray:
    import dask.array as da

    log.info(f"Averaging intensities over {len(cells)} cells")

    tree = shapely.STRtree(cells)

    intensities = np.zeros((len(cells), len(xarr.coords["c"])))
    areas = np.zeros(len(cells))

    def func(x, block_info=None):
        if block_info is not None:
            (ymin, ymax), (xmin, xmax) = block_info[0]["array-location"][1:]
            patch = box(xmin, ymin, xmax, ymax)
            intersections = tree.query(patch, predicate="intersects")

            for index in intersections:
                poly = cells[index]
                bounds = pixel_outer_bounds(poly.bounds)

                sub_image = x[
                    :,
                    max(bounds[1] - ymin, 0) : bounds[3] - ymin,
                    max(bounds[0] - xmin, 0) : bounds[2] - xmin,
                ]

                if sub_image.shape[1] == 0 or sub_image.shape[2] == 0:
                    continue

                mask = rasterize(poly, sub_image.shape[1:], bounds)

                intensities[index] += np.sum(sub_image * mask, axis=(1, 2))
                areas[index] += np.sum(mask)
        return da.zeros_like(x)

    xarr.data.rechunk({0: -1}).map_blocks(func).compute()

    return intensities / areas[:, None].clip(1)
