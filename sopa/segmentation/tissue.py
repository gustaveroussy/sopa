import logging

import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel

from .._constants import ROI
from .._sdata import get_key

log = logging.getLogger(__name__)


def hsv_otsu(
    sdata: SpatialData,
    image_key: str | None = None,
    level: int = -1,
    blur_k: int = 5,
    open_k: int = 5,
    close_k: int = 5,
    drop_threshold: int = 0.01,
):
    """Perform WSI tissue segmentation. The resulting ROI are saved as shapes.

    Args:
        sdata: A `SpatialData` object representing an H&E image
        image_key: Optional key of the H&E image
        level: Level of the multiscale image on which
        blur_k: TODO@stergios
        open_k: TODO@stergios
        close_k: TODO@stergios
        drop_threshold: TODO@stergios
    """
    import cv2

    assert (
        ROI.KEY not in sdata.shapes
    ), f"sdata['{ROI.KEY}'] was already existing, but tissue segmentation is run on top. Delete the shape(s) first."

    image_key = get_key(sdata, "images", image_key)

    level_keys = list(sdata[image_key].keys())
    image: xr.DataArray = next(iter(sdata[image_key][level_keys[level]].values()))

    thumbnail = np.array(image.transpose("y", "x", "c"))
    thumbnail_hsv = cv2.cvtColor(thumbnail, cv2.COLOR_RGB2HSV)
    thumbnail_hsv_blurred = cv2.medianBlur(thumbnail_hsv[:, :, 1], blur_k)
    _, mask = cv2.threshold(thumbnail_hsv_blurred, 0, 255, cv2.THRESH_OTSU + cv2.THRESH_BINARY)

    mask_open = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((open_k, open_k), np.uint8))
    mask_open_close = cv2.morphologyEx(
        mask_open, cv2.MORPH_CLOSE, np.ones((close_k, close_k), np.uint8)
    )

    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask_open_close, 4, cv2.CV_32S)

    contours = []
    for i in range(1, num_labels):
        if stats[i, 4] > drop_threshold * np.prod(mask_open_close.shape):
            cc = cv2.findContours(
                np.array(labels == i, dtype="uint8"),
                cv2.RETR_TREE,
                cv2.CHAIN_APPROX_NONE,
            )[0][0]
            c_closed = np.array(list(cc) + [cc[0]])
            contours.extend([c_closed.squeeze()])

    geo_df = gpd.GeoDataFrame(geometry=[Polygon(contour) for contour in contours])
    geo_df = ShapesModel.parse(
        geo_df,
        transformations=image.attrs["transform"],
    )

    sdata.add_shapes(ROI.KEY, geo_df)

    log.info(f"Tissue segmentation saved in sdata['{ROI.KEY}']")
