from __future__ import annotations

import logging
import warnings

import geopandas as gpd
import numpy as np
from datatree import DataTree
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from xarray import DataArray

from .._constants import SopaAttrs, SopaKeys
from ..utils import add_spatial_element, get_spatial_element
from .shapes import to_valid_polygons

log = logging.getLogger(__name__)

AVAILABLE_MODES = ["hsv_otsu", "staining"]


def tissue_segmentation(
    sdata: SpatialData,
    image_key: str | None = None,
    level: int = -1,
    blur_k: int = 5,
    open_k: int = 5,
    close_k: int = 5,
    drop_threshold: int = 0.01,
    mode: str | None = None,
    key_added: str = SopaKeys.ROI,
):
    """Perform WSI tissue segmentation. The resulting regions-of-interests (ROI) are saved as shapes.

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
        level: Level of the multiscale image on which the segmentation will be performed (if the image is a `DataTree`)
        blur_k: The kernel size of the median bluring operation
        open_k: The kernel size of the morphological openning operation
        close_k: The kernel size of the morphological closing operation
        drop_threshold: Segments that cover less area than `drop_threshold`*100% of the number of pixels of the image will be removed
        mode: Two modes are available: `hsv_otsu` (for H&E data) and `staining`. By default, `hsv_otsu` is used only if there are exactly 3 channels.
        key_added: Name of the spatial element that will be added, containing the segmented tissue polygons.
    """
    assert mode is None or mode in AVAILABLE_MODES, f"`mode` argument should be one of {AVAILABLE_MODES}"

    if key_added in sdata.shapes:
        log.warning(f"sdata['{key_added}'] was already existing, but tissue segmentation is run on top")

    image_key, image = get_spatial_element(
        sdata.images,
        key=image_key or sdata.attrs.get(SopaAttrs.TISSUE_SEGMENTATION),
        return_key=True,
    )

    if isinstance(image, DataTree):
        level_keys = list(image.keys())
        image: DataArray = next(iter(image[level_keys[level]].values()))

    geo_df = TissueSegmentation(image, blur_k, open_k, close_k, drop_threshold).get_polygons(mode)

    if not len(geo_df):
        log.warning(
            "No polygon has been found after tissue segmentation. "
            "Check that there is some tissue in the image, or consider updating the function parameters."
        )
        return

    geo_df = to_valid_polygons(geo_df)
    geo_df = ShapesModel.parse(geo_df, transformations=get_transformation(image, get_all=True).copy())

    add_spatial_element(sdata, key_added, geo_df)


class TissueSegmentation:
    def __init__(self, image: DataArray, blur_k: int, open_k: int, close_k: int, drop_threshold: int):
        self.image = image
        self.blur_k = blur_k
        self.open_k = open_k
        self.close_k = close_k
        self.drop_threshold = drop_threshold

    def _get_default_mode(self) -> str:
        return "hsv_otsu" if self.image.sizes["c"] == 3 else "staining"

    def get_polygons(self, mode: str | None) -> gpd.GeoDataFrame:
        return getattr(self, mode or self._get_default_mode())()

    def hsv_otsu(self) -> gpd.GeoDataFrame:
        import cv2

        thumbnail = np.array(self.image.transpose("y", "x", "c"))

        assert thumbnail.shape[2] == 3, "The image should be in RGB color space"

        if thumbnail.shape[0] * thumbnail.shape[1] > 1e8:
            log.warning(
                "Tissue segmentation is computationally expensive for large images. Consider using a smaller image, or set the `level` parameter."
            )

        thumbnail_hsv = cv2.cvtColor(thumbnail, cv2.COLOR_RGB2HSV)
        thumbnail_hsv_blurred = cv2.medianBlur(thumbnail_hsv[:, :, 1], self.blur_k)
        _, mask = cv2.threshold(thumbnail_hsv_blurred, 0, 255, cv2.THRESH_OTSU + cv2.THRESH_BINARY)

        mask_open = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((self.open_k, self.open_k), np.uint8))
        mask_open_close = cv2.morphologyEx(mask_open, cv2.MORPH_CLOSE, np.ones((self.close_k, self.close_k), np.uint8))

        num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask_open_close, 4, cv2.CV_32S)

        contours = []
        for i in range(1, num_labels):
            if stats[i, 4] > self.drop_threshold * np.prod(mask_open_close.shape):
                cc = cv2.findContours(
                    np.array(labels == i, dtype="uint8"),
                    cv2.RETR_TREE,
                    cv2.CHAIN_APPROX_NONE,
                )[
                    0
                ][0]
                c_closed = np.array(list(cc) + [cc[0]])
                contours.extend([c_closed.squeeze()])

        return gpd.GeoDataFrame(geometry=[Polygon(contour) for contour in contours])

    def staining(self) -> gpd.GeoDataFrame:
        raise NotImplementedError


def hsv_otsu(*args, **kwargs):
    warnings.warn(
        "This function is deprecated and will be removed in late 2024. Use `sopa.tissue_segmentation` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    tissue_segmentation(*args, **kwargs)
