import logging
from math import ceil, floor

import numpy as np
import shapely
import shapely.affinity
from shapely.geometry import MultiPolygon, Polygon
from tqdm import tqdm

log = logging.getLogger(__name__)


def solve_conflicts(
    cells: list[Polygon],
    threshold: float = 0.5,
    patch_indices: np.ndarray | None = None,
    return_indices: bool = False,
) -> np.ndarray[Polygon]:
    """Resolve segmentation conflicts (i.e. overlap) after running segmentation on patches

    Args:
        cells: List of cell polygons
        threshold: When two cells are overlapping, we look at the area of intersection over the area of the smallest cell. If this value is higher than the `threshold`, the cells are merged
        patch_indices: Patch from which each cell belongs.
        return_indices: If `True`, returns also the cells indices.

    Returns:
        Array of resolved cells polygons
    """
    cells = list(cells)
    n_cells = len(cells)
    resolved_indices = np.arange(n_cells)

    assert n_cells > 0, f"No cells was segmented, cannot continue"

    tree = shapely.STRtree(cells)
    conflicts = tree.query(cells, predicate="intersects")

    if patch_indices is not None:
        conflicts = conflicts[:, patch_indices[conflicts[0]] != patch_indices[conflicts[1]]].T
    else:
        conflicts = conflicts[:, conflicts[0] != conflicts[1]].T

    for i1, i2 in tqdm(conflicts, desc="Resolving conflicts"):
        resolved_i1: int = resolved_indices[i1]
        resolved_i2: int = resolved_indices[i2]
        cell1, cell2 = cells[resolved_i1], cells[resolved_i2]

        intersection = cell1.intersection(cell2).area
        if intersection >= threshold * min(cell1.area, cell2.area):
            resolved_indices[np.isin(resolved_indices, [resolved_i1, resolved_i2])] = len(cells)
            cell: Polygon = cell1.union(cell2).buffer(0)
            if cell.interiors:
                cell = Polygon(list(cell.exterior.coords))

            cells.append(cell)

    unique_indices = np.unique(resolved_indices)
    unique_cells = np.array(cells)[unique_indices]

    if return_indices:
        return unique_cells, np.where(unique_indices < n_cells, unique_indices, -1)

    return unique_cells


def _find_contours(cell_mask: np.ndarray) -> MultiPolygon:
    import cv2

    contours, _ = cv2.findContours(cell_mask, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    return MultiPolygon(
        [Polygon(contour[:, 0, :]) for contour in contours if contour.shape[0] >= 4]
    )


def _ensure_polygon(cell: Polygon | MultiPolygon) -> Polygon:
    if isinstance(cell, Polygon):
        return cell

    if isinstance(cell, MultiPolygon):
        return max(cell.geoms, key=lambda polygon: polygon.area)

    log.warn(f"Invalid Polygon type: {type(cell)}. It will not be kept")
    return None


def _smoothen_cell(cell: MultiPolygon, smooth_radius: float, tolerance: float) -> Polygon | None:
    """Smoothen a cell polygon

    Args:
        cell_id: ID of the cell to geometrize
        smooth_radius: radius used to smooth the cell polygon
        tolerance: tolerance used to simplify the cell polygon

    Returns:
        Shapely polygon representing the cell, or `None` if the cell was empty after smoothing
    """
    cell = MultiPolygon([geom for geom in cell.geoms if not geom.buffer(-smooth_radius).is_empty])
    cell = cell.buffer(smooth_radius).buffer(-smooth_radius).simplify(tolerance).buffer(0)

    return None if cell.is_empty else _ensure_polygon(cell)


def _default_tolerance(mean_radius: float) -> float:
    if mean_radius < 10:
        return 0.5
    if mean_radius < 20:
        return 1
    return 2


def geometrize(
    mask: np.ndarray, tolerance: float | None = None, smooth_radius_ratio: float = 0.1
) -> list[Polygon]:
    """Convert a cells mask to multiple `shapely` geometries. Inspired from https://github.com/Vizgen/vizgen-postprocessing

    Args:
        mask: A cell mask. Non-null values correspond to cell ids
        tolerance: Tolerance parameter used by `shapely` during simplification. By default, define the tolerance automatically.

    Returns:
        List of `shapely` polygons representing each cell ID of the mask
    """
    cells = [
        _find_contours((mask == cell_id).astype("uint8")) for cell_id in range(1, mask.max() + 1)
    ]

    mean_radius = np.sqrt(np.array([cell.area for cell in cells]) / np.pi).mean()
    smooth_radius = mean_radius * smooth_radius_ratio

    if tolerance is None:
        tolerance = _default_tolerance(mean_radius)

    cells = [_smoothen_cell(cell, smooth_radius, tolerance) for cell in cells]
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

    coords = np.array(cell_translated.exterior.coords)[None, :].astype(np.int32)
    return cv2.fillPoly(np.zeros((ymax - ymin, xmax - xmin), dtype=np.int8), coords, color=1)
