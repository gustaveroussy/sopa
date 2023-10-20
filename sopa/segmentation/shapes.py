import logging
from math import ceil, floor

import numpy as np
import pandas as pd
import shapely
import shapely.affinity
from shapely.geometry import MultiPolygon, Point, Polygon

log = logging.getLogger(__name__)


def solve_conflicts(
    cells: list[Polygon],
    threshold: float = 0.5,
    patch_indices: np.ndarray | None = None,
    return_indices: bool = False,
) -> np.ndarray[Polygon]:
    n_cells = len(cells)
    resolved_indices = np.arange(n_cells)

    tree = shapely.STRtree(cells)
    conflicts = tree.query(cells, predicate="intersects")

    if patch_indices is not None:
        conflicts = conflicts[:, patch_indices[conflicts[0]] != patch_indices[conflicts[1]]].T
    else:
        conflicts = conflicts[:, conflicts[0] != conflicts[1]].T

    for i1, i2 in conflicts:
        resolved_i1, resolved_i2 = resolved_indices[i1], resolved_indices[i2]
        cell1, cell2 = cells[resolved_i1], cells[resolved_i2]

        intersection = cell1.intersection(cell2).area
        if intersection / min(cell1.area, cell2.area) >= threshold:
            resolved_indices[np.isin(resolved_indices, [resolved_i1, resolved_i2])] = len(cells)
            cells.append(cell1.union(cell2).buffer(0))

    unique_indices = np.unique(resolved_indices)
    unique_cells = np.array(cells)[unique_indices]

    if return_indices:
        return unique_cells, np.where(unique_indices < n_cells, unique_indices, -1)

    return unique_cells


def expand(cells: list[Polygon], expand_radius: float) -> list[Polygon]:
    return [cell.buffer(expand_radius) for cell in cells]


def smooth(cell: Polygon, dist: float, tolerance: float) -> Polygon:
    return cell.buffer(-dist).buffer(-2 * dist).buffer(-dist).simplify(tolerance)


def _find_contours(cell_mask: np.ndarray) -> list[Polygon]:
    import cv2

    contours, _ = cv2.findContours(cell_mask, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    return [Polygon(contour[:, 0, :]) for contour in contours if contour.shape[0] >= 4]


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
    cell = MultiPolygon(polygons)
    cell = cell.buffer(smooth_radius).buffer(-smooth_radius).simplify(tolerance)

    if cell.is_empty:
        return None

    if isinstance(cell, Polygon):
        return cell

    log.warn(
        f"""Geometry index {cell_id} is composed of {len(cell.geoms)} polygons of areas: {[p.area for p in cell.geoms]}. Only the polygon corresponding to the largest area will be kept"""
    )

    return max(cell.geoms, key=lambda polygon: polygon.area)


def geometrize(mask: np.ndarray, smooth_radius: int = 5, tolerance: float = 2) -> list[Polygon]:
    cells = [
        _geometrize_cell(mask, cell_id, smooth_radius, tolerance)
        for cell_id in range(1, mask.max() + 1)
    ]
    return [cell for cell in cells if cell is not None]


def pixel_outer_bounds(bounds: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return [floor(bounds[0]), floor(bounds[1]), ceil(bounds[2]) + 1, ceil(bounds[3]) + 1]


def rasterize(
    cell: Polygon, shape: tuple[int, int], xy_min: tuple[int, int] = [0, 0]
) -> np.ndarray:
    """Transform a cell polygon into a numpy array with value 1 where the polygon touches a pixel, else 0.

    Args:
        cell: Cell polygon to rasterize.
        shape: Image shape as a tuple (y, x).
        xy_min: Tuple containing the origin of the image [x0, y0].

    Returns:
        The mask array.
    """
    import cv2

    xmin, ymin, xmax, ymax = [xy_min[0], xy_min[1], xy_min[0] + shape[1], xy_min[1] + shape[0]]

    cell_translated = shapely.affinity.translate(cell, -xmin, -ymin)
    coords = np.array(cell_translated.boundary.coords)[None, :].astype(np.int32)
    return cv2.fillPoly(np.zeros((ymax - ymin, xmax - xmin), dtype=np.int8), coords, color=1)


def where_transcripts_inside_patch(patch: Polygon, partition: pd.DataFrame) -> np.ndarray:
    points = partition[["x", "y"]].apply(Point, axis=1)
    tree = shapely.STRtree(points)
    indices = tree.query(patch, predicate="intersects")
    where = np.full(len(partition), False)
    where[indices] = True
    return where
