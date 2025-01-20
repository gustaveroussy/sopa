import logging
from enum import Enum

import geopandas as gpd
import numpy as np
from shapely.geometry import GeometryCollection, Polygon, box
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from xarray import DataArray, DataTree

from .._constants import SopaAttrs, SopaKeys
from ..utils import add_spatial_element, get_spatial_element
from .shapes import expand_radius, to_valid_polygons

log = logging.getLogger(__name__)


def tissue(
    sdata: SpatialData,
    image_key: str | None = None,
    level: int = -1,
    mode: str | None = None,
    expand_radius_ratio: float = 0.05,
    channel: str | None = None,
    clip_parameters: tuple[float, float] = (0.9, 5),
    blur_kernel_size: int = 5,
    open_kernel_size: int = 5,
    close_kernel_size: int = 5,
    drop_threshold: float = 0.01,
    key_added: str = SopaKeys.ROI,
):
    """Perform a contouring of the tissue (i.e., "tissue-segmentation"). The resulting regions-of-interest(s) are saved as shapes in the `SpatialData` object. There are two
    modes available: `saturation` and `staining`. The `saturation` mode is used for H&E data, while the `staining` mode
    is used on staining images (more details below).

    !!! info "Saturation mode"
        This segmentation method first transforms the image from RBG color space to HSV. Then,
        on the basis of the saturation channel, a median blurring is applied with an element of size `blur_kernel_size`
        before running the Otsu method. Then a morphological opening and closing are applied as a prostprocessing
        step with square elements of size `open_kernel_size` and `close_kernel_size`. Lastly, the connected components
        with size less than `drop_threshold * number_of_pixel_of_the_image` are removed, and the
        rest are converted into polygons.

    !!! info "Staining mode"
        Instead of extracting the saturation channel, the image is converted to a grayscale image by taking the maximum
        value of all channels (or the specified channel, if `"channel"` is given). The rest of the steps are the same as in the saturation mode.

    Args:
        sdata: A `SpatialData` object representing an H&E image
        image_key: Optional key of the H&E image
        level: Level of the multiscale image on which the segmentation will be performed (if the image is a `DataTree`)
        mode: Two modes are available: `saturation` (for H&E data) and `staining`. By default, `saturation` is used only if there are exactly 3 channels.
        expand_radius_ratio: The ratio of the radius of the polygons that will be expanded.
        channel: The channel to use for the `staining` mode. If `None`, the maximum value of all channels is used.
        clip_parameters: Parameters used to get the threshold used to clip the image before converting it to an 8-bit image (only used in "staining" mode). The first parameter is the quantile, and the second is the divisor. By default, the threshold is the 90th quantile divided by 5.
        blur_kernel_size: The kernel size of the median bluring operation
        open_kernel_size: The kernel size of the morphological openning operation
        close_kernel_size: The kernel size of the morphological closing operation
        drop_threshold: Segments that cover less area than a ratio of `drop_threshold` of the number of pixels of the image will be removed
        key_added: Name of the spatial element that will be added, containing the segmented tissue polygons.
    """
    image, mode = _get_image_and_mode(sdata, image_key, mode, channel)

    if key_added in sdata.shapes:
        log.warning(f"sdata['{key_added}'] was already existing, but tissue segmentation is run on top")

    if isinstance(image, DataTree):
        level_keys = list(image.keys())
        image: DataArray = next(iter(image[level_keys[level]].values()))

    geo_df = TissueSegmentation(
        image=image,
        blur_kernel_size=blur_kernel_size,
        open_kernel_size=open_kernel_size,
        close_kernel_size=close_kernel_size,
        drop_threshold=drop_threshold,
        channel=channel,
        clip_parameters=clip_parameters,
    ).get_polygons(mode)

    if not len(geo_df):
        log.warning(
            "No polygon has been found after tissue segmentation. "
            "Check that there is some tissue in the image, or consider updating the function parameters."
        )
        return

    geo_df = expand_radius(geo_df, expand_radius_ratio)
    geo_df = to_valid_polygons(geo_df)
    geo_df = ShapesModel.parse(geo_df, transformations=get_transformation(image, get_all=True).copy())

    add_spatial_element(sdata, key_added, geo_df)


class AvailableModes(Enum):
    SATURATION = "saturation"
    STAINING = "staining"

    @classmethod
    def check_available(cls, mode: str | None):
        if mode is not None:
            available_modes = list(map(lambda c: c.value, cls))
            assert mode in available_modes, f"Mode '{mode}' not available. Available modes are {available_modes}"


class TissueSegmentation:
    SIZE_THRESHOLD_LARGE_IMAGE = 1e8

    def __init__(
        self,
        image: DataArray,
        blur_kernel_size: int,
        open_kernel_size: int,
        close_kernel_size: int,
        drop_threshold: float,
        channel: str | None,
        clip_parameters: tuple[float, float],
    ):
        self.image = image
        self.blur_kernel_size = blur_kernel_size
        self.open_kernel_size = open_kernel_size
        self.close_kernel_size = close_kernel_size
        self.drop_threshold = drop_threshold
        self.channel = channel
        self.clip_parameters = clip_parameters

        if image.sizes["y"] * image.sizes["x"] > self.SIZE_THRESHOLD_LARGE_IMAGE:
            log.warning(
                "Tissue segmentation is computationally expensive for large images. "
                "Consider using a smaller image, or set the `level` parameter."
            )

    def get_polygons(self, mode: str) -> gpd.GeoDataFrame:
        thumbnail_2d = getattr(self, mode)()

        return self.otsu(thumbnail_2d)

    def staining(self) -> np.ndarray:
        if self.channel is None:
            thumbnail_2d = np.array(self.image.max(dim="c").transpose("y", "x"))
        else:
            thumbnail_2d = np.array(self.image.sel(c=self.channel).transpose("y", "x"))

        threshold = np.quantile(thumbnail_2d, self.clip_parameters[0]) / self.clip_parameters[1]
        thumbnail_2d = ((thumbnail_2d / threshold).clip(0, 1) * 255).astype(np.uint8)

        return thumbnail_2d

    def saturation(self) -> np.ndarray:
        assert self.image.sizes["c"] == 3, "The image should be in RGB color space"

        import cv2

        thumbnail = np.array(self.image.transpose("y", "x", "c"))

        assert thumbnail.dtype == np.uint8, "In 'saturation' mode, the image should have the uint8 dtype"

        thumbnail = cv2.cvtColor(thumbnail, cv2.COLOR_RGB2HSV)
        return thumbnail[:, :, 1]  # saturation channel

    def otsu(self, thumbnail_2d: np.ndarray) -> gpd.GeoDataFrame:
        """Peform tissue segmentation using Otsu's method.

        Args:
            thumbnail_2d: A low-resolution 2D image (YX dimensions) of the tissue.

        Returns:
            A GeoDataFrame containing the segmented polygon(s).
        """
        import cv2

        thumbnail_blurred = cv2.medianBlur(thumbnail_2d, self.blur_kernel_size)
        _, mask = cv2.threshold(thumbnail_blurred, 0, 255, cv2.THRESH_OTSU + cv2.THRESH_BINARY)

        mask_open = cv2.morphologyEx(
            mask, cv2.MORPH_OPEN, np.ones((self.open_kernel_size, self.open_kernel_size), np.uint8)
        )
        mask_open_close = cv2.morphologyEx(
            mask_open, cv2.MORPH_CLOSE, np.ones((self.close_kernel_size, self.close_kernel_size), np.uint8)
        )

        num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask_open_close, 4, cv2.CV_32S)

        contours = []
        for i in range(1, num_labels):
            if stats[i, 4] > self.drop_threshold * np.prod(mask_open_close.shape):
                contours_, _ = cv2.findContours(
                    np.array(labels == i, dtype=np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE
                )
                closed_contours = np.array(list(contours_[0]) + [contours_[0][0]])
                contours.extend([closed_contours.squeeze()])

        return gpd.GeoDataFrame(geometry=[Polygon(contour) for contour in contours])


def hsv_otsu(
    sdata: SpatialData,
    image_key: str | None = None,
    level: int = -1,
    blur_k: int = 5,
    open_k: int = 5,
    close_k: int = 5,
    drop_threshold: float = 0.01,
):
    log.warning(
        "The hsv_otsu function is deprecated and will be removed in sopa==2.1.0. Use `sopa.segmentation.tissue` instead."
    )
    tissue(
        sdata,
        image_key=image_key,
        mode=AvailableModes.SATURATION.value,
        level=level,
        blur_kernel_size=blur_k,
        open_kernel_size=open_k,
        close_kernel_size=close_k,
        drop_threshold=drop_threshold,
    )


def _get_image_and_mode(
    sdata: SpatialData, image_key: str | None, mode: str | None, channel: str | None
) -> tuple[DataArray | DataTree, str]:
    AvailableModes.check_available(mode)
    assert (
        channel is None or mode != AvailableModes.SATURATION.value
    ), "The `channel` argument is only used in the `staining` mode"

    if channel is not None:
        mode = AvailableModes.STAINING.value

    if mode == AvailableModes.STAINING.value:
        image = get_spatial_element(sdata.images, key=image_key or sdata.attrs.get(SopaAttrs.CELL_SEGMENTATION))
        return image, mode

    if mode == AvailableModes.SATURATION.value:
        image = get_spatial_element(sdata.images, key=image_key or sdata.attrs.get(SopaAttrs.TISSUE_SEGMENTATION))
        return image, mode

    if image_key is not None:
        image = get_spatial_element(sdata.images, key=image_key)
        return image, _get_default_mode(image)

    # now, we have no image_key, no mode, no channel

    image_key, image = get_spatial_element(
        sdata.images,
        key=sdata.attrs.get(SopaAttrs.TISSUE_SEGMENTATION) or sdata.attrs.get(SopaAttrs.CELL_SEGMENTATION),
        return_key=True,
    )
    mode = _get_default_mode(image)

    log.info(f"Using {image_key=} and {mode=} as default")

    return image, mode


def _get_default_mode(image: DataArray | DataTree) -> str:
    n_channels = image.sizes["c"] if isinstance(image, DataArray) else next(iter(image.values())).sizes["c"]
    return AvailableModes.SATURATION.value if n_channels == 3 else AvailableModes.STAINING.value


def shapes_bounding_box(sdata: SpatialData, shape_key: str, key_added: str = SopaKeys.ROI):
    """Compute the bounding box of a shape layer and save it as a new shape in the `SpatialData` object.
    For instance, this can be used on VisiumHD data to use the bounding box of the bins as the ROI.

    Args:
        sdata: A `SpatialData` object
        shape_key: Key of the shape that will be used to compute the bounding box
        key_added: Name of the new shape that will be added to the `SpatialData` object
    """
    shapes = get_spatial_element(sdata.shapes, key=shape_key)

    if key_added in sdata.shapes:
        log.warning(f"sdata['{key_added}'] was already existing, but shapes_bounding_box is run on top")

    bounding_box = GeometryCollection(shapes.geometry.values)
    bounding_box = gpd.GeoDataFrame(geometry=[box(*bounding_box.bounds)])

    geo_df = ShapesModel.parse(bounding_box, transformations=get_transformation(shapes, get_all=True).copy())

    add_spatial_element(sdata, key_added, geo_df)
