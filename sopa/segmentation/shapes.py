import logging
from math import ceil, floor

import geopandas as gpd
import numpy as np
import shapely
import shapely.affinity
from shapely.errors import GEOSException
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon

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


def expand_radius(
    geo_df: gpd.GeoDataFrame, expand_radius_ratio: float | None, no_overlap: bool = False, inplace: bool = False
) -> gpd.GeoDataFrame:
    """Expand the radius of the cells by a given ratio.

    Args:
        geo_df: A GeoDataFrame containing the cells or shapes.
        expand_radius_ratio: Ratio to expand the cells polygons for channels averaging. For instance, a ratio of 0.5 expands the shape radius by 50%. If `None`, doesn't expand cells.
        no_overlap: *Experimental feature*: if `True`, ensures that the expanded cells do not overlap by computing the Voronoi diagram of the centroids of the cells.
        inplace: If `True`, modifies the input GeoDataFrame in place. If `False`, returns a new GeoDataFrame.

    Returns:
        A GeoDataFrame with the expanded cells.
    """
    if not inplace:
        geo_df = geo_df.copy()

    if not expand_radius_ratio:
        return geo_df

    expand_radius_ = expand_radius_ratio * np.mean(np.sqrt(geo_df.area / np.pi))
    geo_df.geometry = geo_df.buffer(expand_radius_)

    if no_overlap:
        log.warning("Computing Voronoi polygons to ensure no overlap between shapes is still experimental.")
        geo_df.geometry = remove_overlap(geo_df, as_gdf=False)

    return geo_df


def remove_overlap(geo_df: gpd.GeoDataFrame, as_gdf: bool = True) -> gpd.GeoSeries | gpd.GeoDataFrame:
    """Remove overlapping areas from a GeoDataFrame by computing the Voronoi polygons of the shapes.

    Args:
        geo_df: A GeoDataFrame containing the shapes.
        as_gdf: Whether to return a GeoDataFrame or a GeoSeries.

    Returns:
        A GeoSeries or GeoDataFrame with the overlapping areas removed.
    """
    geo_df["_index"] = geo_df.index  # to keep track of the index after the overlay

    overlay = geo_df.overlay(geo_df, how="intersection")
    overlap = overlay[overlay["_index_1"] != overlay["_index_2"]].union_all()

    if overlap.is_empty:
        return geo_df.geometry

    shapes_no_overlap = geo_df.difference(overlap).buffer(-1e-5)
    _voronoi = voronoi_frames(shapes_no_overlap)

    geometry = geo_df.intersection(_voronoi)

    if not as_gdf:
        return geometry

    geo_df = geo_df.copy()
    geo_df.geometry = geometry
    return geo_df


def voronoi_frames(
    geometry: gpd.GeoSeries | gpd.GeoDataFrame,
    clip: str | shapely.Geometry | None = "bounding_box",
    shrink: float = 0,
    segment: float = 0,
    grid_size: float = 1e-5,
) -> gpd.GeoSeries:
    """
    Copied and simplified from https://pysal.org/libpysal/generated/libpysal.cg.voronoi_frames.html
    """
    # Check if the input geometry is in a geographic CRS
    if geometry.crs and geometry.crs.is_geographic:
        raise ValueError(
            "Geometry is in a geographic CRS. "
            "Use 'GeoSeries.to_crs()' to re-project geometries to a "
            "projected CRS before using voronoi_polygons.",
        )

    # Set precision of the input geometry (avoids GEOS precision issues)
    objects: gpd.GeoSeries = shapely.set_precision(geometry.geometry.copy(), grid_size)

    geom_types = objects.geom_type
    mask_poly = geom_types.isin(["Polygon", "MultiPolygon"])
    mask_line = objects.geom_type.isin(["LineString", "MultiLineString"])

    if mask_poly.any():
        # Shrink polygons if required
        if shrink != 0:
            objects[mask_poly] = objects[mask_poly].buffer(-shrink, cap_style=2, join_style=2)
        # Segmentize polygons if required
        if segment != 0:
            objects.loc[mask_poly] = shapely.segmentize(objects[mask_poly], segment)

    if mask_line.any():
        if segment != 0:
            objects.loc[mask_line] = shapely.segmentize(objects[mask_line], segment)

        # Remove duplicate coordinates from lines
        objects.loc[mask_line] = (
            objects.loc[mask_line]
            .get_coordinates(index_parts=True)
            .drop_duplicates(keep=False)
            .groupby(level=0)
            .apply(shapely.multipoints)
            .values
        )

    limit = _get_limit(objects, clip)

    # Compute Voronoi polygons
    voronoi = shapely.voronoi_polygons(shapely.GeometryCollection(objects.values), extend_to=limit)
    # Get individual polygons out of the collection
    polygons = gpd.GeoSeries(shapely.make_valid(shapely.get_parts(voronoi)), crs=geometry.crs)

    # temporary fix for libgeos/geos#1062
    if not (polygons.geom_type == "Polygon").all():
        polygons = polygons.explode(ignore_index=True)
        polygons = polygons[polygons.geom_type == "Polygon"]

    ids_objects, ids_polygons = polygons.sindex.query(objects, predicate="intersects")

    if mask_poly.any() or mask_line.any():
        # Dissolve polygons
        polygons = polygons.iloc[ids_polygons].groupby(objects.index.take(ids_objects)).agg(_union_with_fallback)
        if geometry.crs is not None:
            polygons = polygons.set_crs(geometry.crs)
    else:
        polygons = polygons.iloc[ids_polygons].reset_index(drop=True)

    # ensure validity as union can occasionally produce invalid polygons that may
    # break the intersection below
    if not polygons.is_valid.all():
        polygons = polygons.make_valid()

    # Clip polygons if limit is provided
    if limit is not None:
        to_be_clipped = polygons.sindex.query(limit.boundary, "intersects")
        polygons.iloc[to_be_clipped] = polygons.iloc[to_be_clipped].intersection(limit)

    return polygons


def _union_with_fallback(arr):
    """
    Coverage union is finnicky with floating point precision and tends to occasionally
    raise an error from within GEOS we have no control over. It is not a data issue
    typically. Falling back to unary union if that happens.
    """
    try:
        r = shapely.coverage_union_all(arr)
    except GEOSException:
        r = shapely.union_all(arr)

    return r


def _get_limit(points, clip: shapely.Geometry | str | None | bool) -> shapely.Geometry | None:
    if isinstance(clip, shapely.Geometry):
        return clip
    if clip is None or clip is False:
        return None
    if clip.lower() == "bounding_box":
        return shapely.box(*points.total_bounds)
    raise ValueError(
        f"Clip type '{clip}' not understood. Try one of the supported options: [None, 'bounding_box', shapely.Polygon]."
    )
