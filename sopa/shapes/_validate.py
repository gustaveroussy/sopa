import logging

import geopandas as gpd
import shapely
from shapely.errors import GEOSException
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon

log = logging.getLogger(__name__)


def to_valid_polygons(geo_df: gpd.GeoDataFrame, simple_polygon: bool = True) -> gpd.GeoDataFrame:
    geo_df.geometry = geo_df.geometry.map(lambda cell: ensure_polygon(cell, simple_polygon))
    return geo_df[~geo_df.is_empty]


def ensure_polygon(
    cell: Polygon | MultiPolygon | GeometryCollection, simple_polygon: bool = True
) -> Polygon | MultiPolygon:
    """Ensures that the provided cell becomes a Polygon

    Args:
        cell: A shapely Polygon or MultiPolygon or GeometryCollection
        simple_polygon: If True, will return a Polygon without holes. Else, allow holes and MultiPolygon.

    Returns:
        The shape as a Polygon, or an empty Polygon if the cell was invalid
    """
    cell = shapely.make_valid(cell)

    if isinstance(cell, Polygon):
        if simple_polygon and cell.interiors:
            cell = Polygon(cell.exterior)
        return cell

    if isinstance(cell, MultiPolygon):
        return max(cell.geoms, key=lambda polygon: polygon.area) if simple_polygon else cell

    if isinstance(cell, GeometryCollection):
        geoms = [geom for geom in cell.geoms if isinstance(geom, Polygon)]

        if len(geoms) > 1 and not simple_polygon:
            return MultiPolygon(geoms)

        if geoms:
            return max(geoms, key=lambda polygon: polygon.area)

        log.warning(f"Found a cell of type {type(cell)} which contains no Polygon geometry")
        return Polygon()

    log.warning(f"Found a cell of unknown type {type(cell)}")
    return Polygon()


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
        cell = cell.buffer(-smooth_radius).buffer(2 * smooth_radius).buffer(-smooth_radius)
        cell = cell.simplify(tolerance)
    except GEOSException:
        log.warning(f"Failed to smoothen cell with {smooth_radius=} and tolerance {tolerance=}.")
        return Polygon()

    return ensure_polygon(cell)


def _default_tolerance(mean_radius: float) -> float:
    if mean_radius < 10:
        return 0.4
    if mean_radius < 20:
        return 1
    return 2
