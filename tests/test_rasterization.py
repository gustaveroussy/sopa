import numpy as np
from shapely.geometry import Polygon

from sopa import shapes

cell = Polygon([[10, 10], [30, 40], [90, 50], [100, 20]])
image_shape = (70, 120)

cell_rectangle = Polygon([[1, 1], [1, 3], [2, 3], [2, 1]])


def test_rasterize():
    mask = shapes.rasterize(cell_rectangle, (5, 5))

    expected = np.array(
        [[0, 0, 0, 0, 0], [0, 1, 1, 0, 0], [0, 1, 1, 0, 0], [0, 1, 1, 0, 0], [0, 0, 0, 0, 0]],
        dtype=np.int8,
    )

    assert (mask == expected).all()


def test_raster_and_vectorize():
    mask = shapes.rasterize(cell, image_shape)

    new_cell = shapes.vectorize(mask, 0, 0).geometry[0]
    new_mask = shapes.rasterize(new_cell, image_shape)

    iou = (mask & new_mask).sum() / (mask | new_mask).sum()

    assert iou > 0.99, "Applying vectorize and then rasterize shouldn't change the mask"


def test_rasterize_cropped():
    mask = shapes.rasterize(cell_rectangle, (2, 4), (1, 3))

    assert mask.sum() == 2
