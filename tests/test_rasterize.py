import numpy as np
from shapely.geometry import Polygon

from sopa.segmentation import shapes

polygon = Polygon([[1, 1], [3, 4], [9, 5], [10, 2]])
image_shape = (7, 12)

rectangle = Polygon([[1, 1], [1, 3], [2, 3], [2, 1]])


def test_rasterize():
    mask = shapes.rasterize(rectangle, (5, 5))

    expected = np.array(
        [[0, 0, 0, 0, 0], [0, 1, 1, 0, 0], [0, 1, 1, 0, 0], [0, 1, 1, 0, 0], [0, 0, 0, 0, 0]],
        dtype=np.int8,
    )

    assert (mask == expected).all()


def test_raster_and_geometrize():
    mask = shapes.rasterize(polygon, image_shape)

    new_poly = shapes.geometrize(mask, 0, 0)[0]
    new_mask = shapes.rasterize(new_poly, image_shape)

    assert (
        mask == new_mask
    ).all(), f"Applying geometrize and then rasterize shouldn't change the mask"


def test_rasterize_cropped():
    mask = shapes.rasterize(rectangle, (2, 4), (1, 3))

    assert mask.sum() == 2
