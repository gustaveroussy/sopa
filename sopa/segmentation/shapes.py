import logging
from math import ceil, floor

import numpy as np
import rasterio.features
import shapely
import shapely.affinity
import xarray as xr
from shapely.geometry import MultiPolygon, Polygon
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


def smooth(poly: Polygon, smooth_radius: int = 5) -> Polygon:
    # Copied from https://github.com/Vizgen/vizgen-postprocessing
    smooth = poly.buffer(-smooth_radius).buffer(smooth_radius * 2).buffer(-smooth_radius)
    return poly if isinstance(smooth, MultiPolygon) else smooth


def extract_polygons(mask: np.ndarray) -> list[Polygon]:
    # Copied from https://github.com/Vizgen/vizgen-postprocessing
    # TODO: do not rely on cv2?
    import cv2

    polys = []

    for cell_id in range(1, mask.max() + 1):
        mask_id = (mask == cell_id).astype("uint8")
        contours, _ = cv2.findContours(mask_id, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

        polys_ = [smooth(Polygon(c[:, 0, :])) for c in contours if c.shape[0] >= 4]
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


def to_chunk_mask(poly: Polygon, bounds: tuple[int, int, int, int]) -> np.ndarray:
    xmin, ymin, xmax, ymax = bounds

    new_poly = shapely.affinity.translate(poly, -xmin, -ymin)
    return rasterio.features.rasterize([new_poly], out_shape=(ymax - ymin, xmax - xmin))


def average_polygon(xarr: xr.DataArray, poly: Polygon) -> np.ndarray:
    bounds = outer_bounds(poly.bounds)

    sub_image = xarr.sel(
        x=slice(bounds[0], bounds[2]), y=slice(bounds[1], bounds[3])
    ).data.compute()

    bounds = update_bounds(bounds, sub_image.shape[1:])
    mask = to_chunk_mask(poly, bounds)
    return np.sum(sub_image * mask, axis=(1, 2)) / np.sum(mask)


def average(xarr: xr.DataArray, polygons: list[Polygon]) -> np.ndarray:
    log.info(f"Averaging intensities over {len(polygons)} polygons")
    return np.stack([average_polygon(xarr, poly) for poly in tqdm(polygons)])
