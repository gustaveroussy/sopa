import logging
from math import ceil, floor

import geopandas as gpd
import numpy as np
import shapely
import shapely.affinity
import skimage
from geopandas import GeoDataFrame
from shapely.geometry import GeometryCollection, LinearRing, MultiPolygon, Polygon
from skimage.draw import polygon
from skimage.measure._regionprops import RegionProperties

log = logging.getLogger(__name__)


def _ensure_polygon(
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
            cell = Polygon(list(cell.exterior.coords))
        return cell

    if isinstance(cell, MultiPolygon):
        return max(cell.geoms, key=lambda polygon: polygon.area) if simple_polygon else cell

    if isinstance(cell, GeometryCollection):
        geoms = [geom for geom in cell.geoms if isinstance(geom, Polygon)]

        if len(geoms) > 1 and not simple_polygon:
            return MultiPolygon(geoms)

        if geoms:
            return max(geoms, key=lambda polygon: polygon.area)

        geoms = [geom for geom in cell.geoms if isinstance(geom, MultiPolygon)]
        geoms = [polygon for multi_polygon in geoms for polygon in multi_polygon.geoms]

        if len(geoms) > 1 and not simple_polygon:
            return MultiPolygon(geoms)

        if geoms:
            return max(geoms, key=lambda polygon: polygon.area)

        log.warning(f"Removing cell of type {type(cell)} as it contains no Polygon geometry")
        return Polygon()

    log.warning(f"Removing cell of unknown type {type(cell)}")
    return Polygon()


def to_valid_polygons(geo_df: gpd.GeoDataFrame, simple_polygon: bool = True) -> gpd.GeoDataFrame:
    geo_df.geometry = geo_df.geometry.map(lambda cell: _ensure_polygon(cell, simple_polygon))
    return geo_df[~geo_df.is_empty]


def _smoothen_cell(cell: MultiPolygon, smooth_radius: float, tolerance: float) -> Polygon:
    """Smoothen a cell polygon

    Args:
        cell_id: MultiPolygon representing a cell
        smooth_radius: radius used to smooth the cell polygon
        tolerance: tolerance used to simplify the cell polygon

    Returns:
        Shapely polygon representing the cell, or an empty Polygon if the cell was empty after smoothing
    """
    cell = cell.buffer(-smooth_radius).buffer(2 * smooth_radius).buffer(-smooth_radius)
    cell = cell.simplify(tolerance)

    return _ensure_polygon(cell)


def _default_tolerance(mean_radius: float) -> float:
    if mean_radius < 10:
        return 0.4
    if mean_radius < 20:
        return 1
    return 2


def geometrize(*args, **kwargs):
    """Alias for `vectorize`"""
    return vectorize(*args, **kwargs)


def vectorize(mask: np.ndarray, tolerance: float | None = None, smooth_radius_ratio: float = 0.1) -> gpd.GeoDataFrame:
    """Convert a cells mask to multiple `shapely` geometries. Inspired from https://github.com/Vizgen/vizgen-postprocessing

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


def _region_props_to_multipolygon(region_props: RegionProperties, allow_holes: bool) -> MultiPolygon:
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

    multipolygon = MultiPolygon([_to_polygon(exterior) for exterior in exteriors])

    yoff, xoff, *_ = region_props.bbox
    return shapely.affinity.translate(multipolygon, xoff - 1, yoff - 1)  # remove padding offset


def _vectorize_mask(mask: np.ndarray, allow_holes: bool = False) -> GeoDataFrame:
    if mask.max() == 0:
        return GeoDataFrame(geometry=[])

    regions = skimage.measure.regionprops(mask)

    return GeoDataFrame(geometry=[_region_props_to_multipolygon(region, allow_holes) for region in regions])


def pixel_outer_bounds(bounds: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return [floor(bounds[0]), floor(bounds[1]), ceil(bounds[2]) + 1, ceil(bounds[3]) + 1]


def rasterize(cell: Polygon | MultiPolygon, shape: tuple[int, int], xy_min: tuple[int, int] = (0, 0)) -> np.ndarray:
    """Transform a cell polygon into a numpy array with value 1 where the polygon touches a pixel, else 0.

    Args:
        cell: Cell polygon to rasterize.
        shape: Image shape as a tuple (y, x).
        xy_min: Tuple containing the origin of the image [x0, y0].

    Returns:
        The mask array.
    """
    xmin, ymin, xmax, ymax = [xy_min[0], xy_min[1], xy_min[0] + shape[1], xy_min[1] + shape[0]]

    cell_translated = shapely.affinity.translate(cell, -xmin, -ymin)
    geoms = cell_translated.geoms if isinstance(cell_translated, MultiPolygon) else [cell_translated]

    shape = (ymax - ymin, xmax - xmin)

    rasterized_image = np.zeros(shape, dtype=np.int8)

    for geom in geoms:
        x, y = geom.exterior.coords.xy
        rr, cc = polygon(y, x, shape)
        rasterized_image[rr, cc] = 1

    return rasterized_image


def expand_radius(geo_df: gpd.GeoDataFrame, expand_radius_ratio: float | None) -> gpd.GeoDataFrame:
    """Expand the radius of the cells by a given ratio.

    Args:
        geo_df: A GeoDataFrame containing the cells or shapes.
        expand_radius_ratio: Ratio to expand the cells polygons for channels averaging. For instance, a ratio of 0.5 expands the shape radius by 50%. If `None`, doesn't expand cells.

    Returns:
        A GeoDataFrame with the expanded cells.
    """
    if not expand_radius_ratio:
        return geo_df

    expand_radius_ = expand_radius_ratio * np.mean(np.sqrt(geo_df.area / np.pi))
    geo_df.geometry = geo_df.buffer(expand_radius_)
    return geo_df
