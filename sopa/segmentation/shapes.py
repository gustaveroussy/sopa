import logging
from math import ceil, floor

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
    import cv2

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


def outer_bounds(bounds: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return [floor(bounds[0]), floor(bounds[1]), ceil(bounds[2]) + 1, ceil(bounds[3]) + 1]


def update_bounds(
    bounds: tuple[int, int, int, int], shape: tuple[int, int]
) -> tuple[int, int, int, int]:
    """Update bound's width and heigh to fit the image

    Args:
        bounds: original bounds (xmin, ymin, xmax, ymax)
        shape: image shapes (dim_y, dim_x)

    Returns:
        Updated bounds
    """
    return (bounds[0], bounds[1], bounds[0] + shape[1], bounds[1] + shape[0])


def rasterize(poly: Polygon, bounds: tuple[int, int, int, int]) -> np.ndarray:
    import cv2

    xmin, ymin, xmax, ymax = bounds

    new_poly = shapely.affinity.translate(poly, -xmin, -ymin)
    coords = np.array(new_poly.boundary.coords)[None, :].astype(np.int32)
    return cv2.fillPoly(np.zeros((ymax - ymin, xmax - xmin), dtype=np.int8), coords, color=1)


def average_polygon(xarr: xr.DataArray, poly: Polygon) -> np.ndarray:
    bounds = outer_bounds(poly.bounds)

    sub_image = xarr.sel(
        x=slice(bounds[0], bounds[2]), y=slice(bounds[1], bounds[3])
    ).data.compute()

    bounds = update_bounds(bounds, sub_image.shape[1:])
    mask = rasterize(poly, bounds)
    return np.sum(sub_image * mask, axis=(1, 2)) / np.sum(mask)


def average(xarr: xr.DataArray, polygons: list[Polygon]) -> np.ndarray:
    log.info(f"Averaging intensities over {len(polygons)} polygons")
    return np.stack([average_polygon(xarr, poly) for poly in tqdm(polygons)])
