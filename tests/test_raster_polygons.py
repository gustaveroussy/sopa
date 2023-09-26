from shapely.geometry import Polygon

from sopa.segmentation import shapes


def test_raster_inverse():
    poly = Polygon([[1, 1], [3, 4], [9, 5], [10, 2]])
    bounds = [0, 0, 12, 7]
    mask = shapes.rasterize(poly, bounds)

    new_poly = shapes.extract_polygons(mask, 0, 0)[0]
    new_mask = shapes.rasterize(new_poly, bounds)

    assert (mask == new_mask).all()
