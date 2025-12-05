import logging
from math import ceil, floor

import geopandas as gpd
import numpy as np
import shapely
import shapely.affinity
import skimage
from geopandas import GeoDataFrame
from shapely.errors import GEOSException
from shapely.geometry import LinearRing, MultiPolygon, Polygon
from skimage.measure._regionprops import RegionProperties

from .. import settings
from .validate import ensure_polygon

log = logging.getLogger(__name__)


def vectorize(mask: np.ndarray, tolerance: float | None = None, smooth_radius_ratio: float = 0.1) -> gpd.GeoDataFrame:
    """Convert a cells mask to multiple `shapely` geometries. Each shape is smoothen and converted to a single polygon. Inspired from [VPT](https://github.com/Vizgen/vizgen-postprocessing).

    Args:
        mask: A cell mask. Non-null values correspond to cell ids
        tolerance: Tolerance parameter used by `shapely` during simplification. By default, define the tolerance automatically.
        smooth_radius_ratio: Ratio of the cell radius used to smooth the cell polygon.

    Returns:
        GeoDataFrame of polygons representing each cell ID of the mask
    """
    max_cells = mask.max()

    if max_cells == 0:
        log.warning("No cell was returned by the segmentation")
        return gpd.GeoDataFrame(geometry=[])

    cells = _vectorize_mask(mask)

    mean_radius = np.sqrt(cells.area / np.pi).mean()
    smooth_radius = mean_radius * smooth_radius_ratio

    tolerance = _default_tolerance(mean_radius) if tolerance is None else tolerance

    cells.geometry = cells.geometry.map(lambda cell: _smoothen_cell(cell, smooth_radius, tolerance))
    cells = cells[~cells.is_empty]

    return cells


def _region_props_to_multipolygon(region_props: RegionProperties, allow_holes: bool) -> Polygon | MultiPolygon:
    mask = np.pad(region_props.image, 1)
    contours = skimage.measure.find_contours(mask, 0.5)

    rings = [LinearRing(contour[:, [1, 0]]) for contour in contours if contour.shape[0] >= 4]

    exteriors = [ring for ring in rings if ring.is_ccw]

    if allow_holes:
        holes = [ring for ring in rings if not ring.is_ccw]

    def _to_polygon(exterior: LinearRing) -> Polygon:
        exterior_poly = Polygon(exterior)

        _holes = None
        if allow_holes:
            _holes = [hole.coords for hole in holes if exterior_poly.contains(Polygon(hole))]

        return Polygon(exterior, holes=_holes)

    polygon = [_to_polygon(exterior) for exterior in exteriors]
    polygon = MultiPolygon(polygon) if len(polygon) > 1 else polygon[0]

    yoff, xoff, *_ = region_props.bbox
    return shapely.affinity.translate(polygon, xoff - 1, yoff - 1)  # remove padding offset


def _vectorize_mask(mask: np.ndarray, allow_holes: bool = False) -> GeoDataFrame:
    if mask.max() == 0:
        return GeoDataFrame(geometry=[])

    regions = skimage.measure.regionprops(mask)

    return GeoDataFrame(geometry=[_region_props_to_multipolygon(region, allow_holes) for region in regions])


def pixel_outer_bounds(bounds: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return [floor(bounds[0]), floor(bounds[1]), ceil(bounds[2]) + 1, ceil(bounds[3]) + 1]


def _smoothen_cell(cell: MultiPolygon, smooth_radius: float, tolerance: float) -> Polygon:
    """Smoothen a cell polygon

    Args:
        cell_id: MultiPolygon representing a cell
        smooth_radius: radius used to smooth the cell polygon
        tolerance: tolerance used to simplify the cell polygon

    Returns:
        Shapely polygon representing the cell, or an empty Polygon if the cell was empty after smoothing
    """
    try:
        if smooth_radius > 0:
            cell = cell.buffer(-smooth_radius).buffer(2 * smooth_radius).buffer(-smooth_radius)
        if tolerance > 0:
            cell = cell.simplify(tolerance)
    except GEOSException:
        log.warning(f"Failed to smoothen cell with {smooth_radius=} and tolerance {tolerance=}.")
        return Polygon()

    return ensure_polygon(cell)


def _default_tolerance(mean_radius: float) -> float:
    return settings.simplification_tolerance or min(mean_radius * 0.03, 1)
