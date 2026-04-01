import logging

import numpy as np
import shapely
import shapely.affinity
from shapely.geometry import MultiPolygon, Polygon
from skimage.draw import polygon

log = logging.getLogger(__name__)


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
