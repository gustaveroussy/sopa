import logging
from math import ceil, floor

import cv2
import numpy as np
import shapely
import shapely.affinity
import xarray as xr
from shapely.geometry import Polygon
from tqdm import tqdm

log = logging.getLogger(__name__)


def solve_conflicts(polygons: list[Polygon], threshold: float = 0.5) -> np.ndarray[Polygon]:
    n_polygons = len(polygons)
    resolved_indices = np.arange(n_polygons)

    tree = shapely.STRtree(polygons)
    conflicts = tree.query(polygons, predicate="intersects")
    conflicts = conflicts[:, conflicts[0] != conflicts[1]].T

    for i1, i2 in conflicts:
        resolved_i1, resolved_i2 = resolved_indices[i1], resolved_indices[i2]
        poly1, poly2 = polygons[resolved_i1], polygons[resolved_i2]

        intersection = poly1.intersection(poly2).area
        if intersection / min(poly1.area, poly2.area) >= threshold:
            resolved_indices[np.isin(resolved_indices, [resolved_i1, resolved_i2])] = len(polygons)
            polygons.append(poly1.union(poly2))

    return np.array(polygons)[np.unique(resolved_indices)]


def expand(polygons: list[Polygon], expand_radius: float) -> list[Polygon]:
    return [polygon.buffer(expand_radius) for polygon in polygons]


def smooth(poly: Polygon, smooth_radius: int, tolerance: float) -> Polygon:
    return poly.buffer(smooth_radius).buffer(-smooth_radius).simplify(tolerance)


def geometrize(mask: np.ndarray, smooth_radius: int = 3, tolerance: float = 2) -> list[Polygon]:
    # Copied from https://github.com/Vizgen/vizgen-postprocessing
    polys = []

    for cell_id in range(1, mask.max() + 1):
        mask_id = (mask == cell_id).astype("uint8")
        contours, _ = cv2.findContours(mask_id, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

        polys_ = [Polygon(c[:, 0, :]) for c in contours if c.shape[0] >= 4]
        polys_ = [smooth(p, smooth_radius, tolerance) for p in polys_]
        polys_ = [p for p in polys_ if not p.is_empty]

        assert len(polys_) <= 1
        polys.extend(polys_)

    return polys


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
    xmin, ymin, xmax, ymax = [xy_min[0], xy_min[1], xy_min[0] + shape[1], xy_min[1] + shape[0]]

    new_poly = shapely.affinity.translate(poly, -xmin, -ymin)
    coords = np.array(new_poly.boundary.coords)[None, :].astype(np.int32)
    return cv2.fillPoly(np.zeros((ymax - ymin, xmax - xmin), dtype=np.int8), coords, color=1)


def average_polygon(xarr: xr.DataArray, poly: Polygon) -> np.ndarray:
    bounds = pixel_outer_bounds(poly.bounds)

    sub_image = xarr.sel(
        x=slice(bounds[0], bounds[2]), y=slice(bounds[1], bounds[3])
    ).data.compute()

    mask = rasterize(poly, sub_image.shape[1:], bounds[:2])
    return np.sum(sub_image * mask, axis=(1, 2)) / np.sum(mask)


def average(xarr: xr.DataArray, polygons: list[Polygon]) -> np.ndarray:
    log.info(f"Averaging intensities over {len(polygons)} polygons")
    return np.stack([average_polygon(xarr, poly) for poly in tqdm(polygons)])
