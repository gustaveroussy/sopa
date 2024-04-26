from __future__ import annotations

import logging

import geopandas as gpd
import numpy as np
import spatialdata
import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel

from .._constants import ROI
from .._sdata import get_intrinsic_cs, get_item, save_shapes

log = logging.getLogger(__name__)


def hsv_otsu(
    sdata: SpatialData,
    image_key: str | None = None,
    level: int = -1,
    blur_k: int = 5,
    open_k: int = 5,
    close_k: int = 5,
    drop_threshold: int = 0.01,
) -> bool:
    """Perform WSI tissue segmentation. The resulting ROIs are saved as shapes.

    !!! info
        This segmentation method first transforms the image from RBG color space to HSV and then
        on the basis of the saturation channel, it performs the rest of the steps.
        As a preprocessing step, a median blurring is applied with an element of size `blur_k`
        before the otsu. Then a morphological opening and closing are applied as a prostprocessing
        step with square elements of size `open_k` and `close_k`. Lastly, the connected components
        with size less than `drop_threshold * number_of_pixel_of_the_image` are removed, and the
        rest are converted into polygons.

    Args:
        sdata: A `SpatialData` object representing an H&E image
        image_key: Optional key of the H&E image
        level: Level of the multiscale image on which the segmentation will be performed
        blur_k: The kernel size of the median bluring operation
        open_k: The kernel size of the morphological openning operation
        close_k: The kernel size of the morphological closing operation
        drop_threshold: Segments that cover less area than `drop_threshold`*100% of the number of pixels of the image will be removed

    Returns:
        `True` if tissue segmentation was successful, else `False` if no polygon was output.
    """
    import cv2

    image_key, image = get_item(sdata, "images", image_key)

    if level == 0 or not isinstance(image, MultiscaleSpatialImage):
        log.warn(
            "Running hsv_otsu on the full image can be slow. We recommend using a MultiscaleSpatialImage"
        )

    if isinstance(image, MultiscaleSpatialImage):
        level_keys = list(image.keys())
        image: xr.DataArray = next(iter(image[level_keys[level]].values()))

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

    polygons = [Polygon(contour) for contour in contours]

    if not len(polygons):
        log.warn(
            "No polygon has been found after tissue segmentation. Check that there is some tissue in the image, or consider updating the segmentation parameters."
        )
        return False

    _save_tissue_segmentation(sdata, polygons, image_key, image)
    return True


def _save_tissue_segmentation(
    sdata: SpatialData, polygons: list[Polygon], image_key: str, image_scale: xr.DataArray
):
    assert (
        ROI.KEY not in sdata.shapes
    ), f"sdata['{ROI.KEY}'] was already existing, but tissue segmentation is run on top. Delete the shape(s) first."

    geo_df = gpd.GeoDataFrame(geometry=polygons)
    geo_df = ShapesModel.parse(
        geo_df,
        transformations=image_scale.attrs["transform"],
    )

    image_cs = get_intrinsic_cs(sdata, sdata[image_key])
    geo_df = spatialdata.transform(
        geo_df, image_scale.attrs["transform"][image_cs], maintain_positioning=True
    )

    sdata.shapes[ROI.KEY] = geo_df
    save_shapes(sdata, ROI.KEY)

    log.info(f"Tissue segmentation saved in sdata['{ROI.KEY}']")
