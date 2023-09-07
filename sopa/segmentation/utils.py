import numpy as np
from shapely.geometry import MultiPolygon, Polygon


def smooth(poly: Polygon, radius: int = 5) -> Polygon:
    smooth = poly.buffer(-radius).buffer(radius * 2).buffer(-radius)
    return poly if isinstance(smooth, MultiPolygon) else smooth


def pad(
    polygon: Polygon, min_vertices: int, max_vertices: int, tolerance: float = 1
) -> np.ndarray:
    n_vertices = len(polygon.exterior.coords)
    assert n_vertices >= min_vertices

    coords = polygon.exterior.coords._coords

    if n_vertices == max_vertices:
        return coords.flatten()

    if n_vertices < max_vertices:
        return np.pad(
            coords, ((0, max_vertices - n_vertices), (0, 0)), mode="edge"
        ).flatten()

    polygon = polygon.simplify(tolerance=tolerance)
    return pad(polygon, min_vertices, max_vertices, tolerance + 1)
