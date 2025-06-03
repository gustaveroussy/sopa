import cv2
import numpy as np
import skimage
from shapely import MultiPolygon
from shapely.geometry import Polygon
from shapely.ops import unary_union
from spatialdata._core.operations.vectorize import _vectorize_mask
from xarray import DataArray, DataTree

import sopa
from sopa.segmentation._tissue import TissueSegmentation, _get_image_and_mode


def compare_cv2_skimage():
    sdata = sopa.io.toy_dataset()

    image, _ = _get_image_and_mode(sdata, None, None, None)

    if isinstance(image, DataTree):
        level_keys = list(image.keys())
        image: DataArray = next(iter(image[level_keys[-1]].values()))

    tissue_seg = TissueSegmentation(
        image=image,
        blur_kernel_size=5,
        open_kernel_size=5,
        close_kernel_size=5,
        drop_threshold=0.01,
        channel=None,
        clip_parameters=(0.9, 5),
    )

    # skimage thumbnail 2D
    thumbnail = np.array(tissue_seg.image.transpose("y", "x", "c"))
    thumbnail = (skimage.color.rgb2hsv(thumbnail) * 255).astype("uint8")
    thumbnail_2d_skimage = thumbnail[:, :, 1]

    # cv2 thumbnail 2D
    thumbnail = np.array(tissue_seg.image.transpose("y", "x", "c"))
    thumbnail = cv2.cvtColor(thumbnail, cv2.COLOR_RGB2HSV)
    thumbnail_2d_cv = thumbnail[:, :, 1]

    # comparison
    assert thumbnail_2d_skimage.max() == thumbnail_2d_cv.max()
    assert ((thumbnail_2d_skimage - thumbnail_2d_cv) ** 2).max() <= 1
    thumbnail_2d = thumbnail_2d_cv

    # skimage thumbnail blurred
    median_footprint = np.ones((tissue_seg.blur_kernel_size, tissue_seg.blur_kernel_size))
    thumbnail_blurred_skimage = skimage.filters.median(thumbnail_2d, footprint=median_footprint)

    # cv2 thumbnail blurred
    thumbnail_blurred_cv = cv2.medianBlur(thumbnail_2d, 5)

    # comparison
    assert (thumbnail_blurred_skimage == thumbnail_blurred_cv).all()
    thumbnail_blurred = thumbnail_blurred_cv

    # skimage labels
    mask = thumbnail_blurred > skimage.filters.threshold_otsu(thumbnail_blurred)

    opening_footprint = np.ones((tissue_seg.open_kernel_size, tissue_seg.open_kernel_size))
    mask_open = skimage.morphology.opening(mask, opening_footprint).astype(np.uint8)

    closing_footprint = np.ones((tissue_seg.close_kernel_size, tissue_seg.close_kernel_size))
    mask_open_close = skimage.morphology.closing(mask_open, closing_footprint).astype(np.uint8)

    labels_skimage = skimage.measure.label(mask_open_close, connectivity=2)

    # cv2 labels
    _, mask = cv2.threshold(thumbnail_blurred, 0, 255, cv2.THRESH_OTSU + cv2.THRESH_BINARY)

    mask_open = cv2.morphologyEx(
        mask, cv2.MORPH_OPEN, np.ones((tissue_seg.open_kernel_size, tissue_seg.open_kernel_size), np.uint8)
    )

    mask_open_close = cv2.morphologyEx(
        mask_open, cv2.MORPH_CLOSE, np.ones((tissue_seg.close_kernel_size, tissue_seg.close_kernel_size), np.uint8)
    )

    num_labels, labels_cv, stats, _ = cv2.connectedComponentsWithStats(mask_open_close, 4, cv2.CV_32S)

    # comparison
    assert (labels_skimage == labels_cv).all()

    # skimage multi polygon
    gdf_ = _vectorize_mask(labels_skimage)
    multi_polygon_skimage = unary_union(MultiPolygon(gdf_.geometry.tolist()))

    # cv2 multi polygon
    contours = []
    for i in range(1, num_labels):
        if stats[i, 4] > tissue_seg.drop_threshold * np.prod(mask_open_close.shape):
            contours_, _ = cv2.findContours(
                np.array(labels_cv == i, dtype=np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE
            )
            closed_contours = np.array([*list(contours_[0]), contours_[0][0]])
            contours.extend([closed_contours.squeeze()])

    multi_polygon_cv = unary_union(MultiPolygon([Polygon(contour) for contour in contours]))

    # comparison
    iou = multi_polygon_skimage.intersection(multi_polygon_cv).area / multi_polygon_skimage.union(multi_polygon_cv).area
    assert 0.97 <= iou <= 1

    print("skimage and cv2 tissue segmentation are relatively equivalent")


if __name__ == "__main__":
    compare_cv2_skimage()
