import logging

import geopandas as gpd
import numpy as np
import shapely
from shapely.errors import GEOSException

from ._validate import ensure_polygon

log = logging.getLogger(__name__)


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
        log.warning(
            "Computing Voronoi polygons to ensure no overlap between shapes is still experimental. It can take 10+ minutes for 100k+ shapes."
        )
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

    del geo_df["_index"]

    if overlap.is_empty:
        return geo_df.geometry

    shapes_no_overlap = geo_df.difference(overlap).buffer(-1e-4)  # to avoid touching polygons on single points
    _voronoi = voronoi_frames(shapes_no_overlap)

    geometry = geo_df.intersection(_voronoi)

    nan_locs = geometry.type.isna()

    if nan_locs.any():
        log.warning(f"Found {nan_locs.sum()} NaN geometries after removing the overlap. Replacing with empty polygons.")
        geometry[nan_locs] = shapely.Polygon()

    geometry = geometry.buffer(1e-4)  # to re-expand the polygons to their original size
    geometry = geometry.map(ensure_polygon)

    if not as_gdf:
        return geometry

    geo_df = geo_df.copy()
    geo_df.geometry = geometry
    return geo_df


def voronoi_frames(
    geometry: gpd.GeoSeries | gpd.GeoDataFrame,
    clip: str | shapely.Geometry | None = "bounding_box",
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
    assert geom_types.isin(["Polygon", "MultiPolygon"]).all(), (
        "Only Polygon and MultiPolygon geometries are supported to remove overlaps."
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

    # Dissolve polygons
    polygons = polygons.iloc[ids_polygons].groupby(objects.index.take(ids_objects)).agg(_union_with_fallback)
    if geometry.crs is not None:
        polygons = polygons.set_crs(geometry.crs)

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
