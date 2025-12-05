import logging

import dask
import numpy as np
from spatialdata import SpatialData
from xarray import DataArray, DataTree

from .. import settings
from ..constants import SopaAttrs, SopaKeys
from ..io.reader._wsi_reader import get_reader
from ..utils import get_spatial_element
from .patches import Patches2D

log = logging.getLogger(__name__)


class TileLoader:
    def __init__(
        self,
        sdata: SpatialData,
        patch_width: int,
        image_key: str | None = None,
        level: int | None = 0,
        magnification: int | None = None,
        patch_overlap: int = 0,
        roi_key: str | None = SopaKeys.ROI,
    ):
        multiscale_image = _get_image_for_inference(sdata, image_key)
        self.image, self.level, self.tile_resize_factor, self.level_downsample = _get_extraction_parameters(
            multiscale_image, level, magnification
        )

        _backend = multiscale_image.attrs.get("backend")
        _path = multiscale_image.attrs.get("path")

        try:
            if settings.native_read_region and _backend is not None:
                self.slide = get_reader(_backend)(_path)
            else:
                self.slide = get_reader("xarray")(multiscale_image)
        except Exception as e:
            log.warning(
                f"Exception raised for '{_backend}' and path '{_path}'. Falling back to xarray reader. Error: {e}"
            )
            self.slide = get_reader("xarray")(multiscale_image)

        self.patch_width = int(patch_width / self.tile_resize_factor)
        self.resized_patch_width = patch_width
        self.patches = Patches2D(sdata, self.image, self.patch_width, patch_overlap=patch_overlap, roi_key=roi_key)

    def _numpy_patch(self, box: tuple[int, int, int, int]) -> np.ndarray:
        """
        Extract a numpy patch from the image given a bounding box
        and pads a patch to a specific width since some patches might be smaller (e.g., on edges)
        """
        image_patch = np.array(
            self.slide.read_region(
                (int(box[0] * self.level_downsample), int(box[1] * self.level_downsample)),
                self.level,
                (box[2] - box[0], box[3] - box[1]),
            )
        )

        pad_x, pad_y = self.patch_width - image_patch.shape[0], self.patch_width - image_patch.shape[1]
        padded_patch = np.pad(image_patch, ((0, pad_x), (0, pad_y), (0, 0)))
        return padded_patch

    def get_batch(self, bboxes: np.ndarray):
        """Retrives a batch of patches using the bboxes"""
        import torch

        delayed_patches = [dask.delayed(self._numpy_patch)(box) for box in bboxes]
        batch = np.array(dask.compute(*delayed_patches))

        batch = torch.tensor(batch, dtype=torch.float32) / 255.0
        batch = batch.movedim(-1, 1)

        if self.tile_resize_factor != 1:
            from torchvision.transforms import Resize

            dim = (self.resized_patch_width, self.resized_patch_width)
            batch = Resize(dim)(batch)

        assert len(batch.shape) == 4

        return batch

    def get_PIL_patch(self, idx: int):
        """Retrives a PIL patch using the index"""
        from PIL import Image

        return Image.fromarray(np.array(self[idx][0].movedim(0, -1) * 255, dtype="uint8"))

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, idx: int | slice):
        bboxes = self.patches.bboxes[idx]
        batch = self.get_batch(np.atleast_2d(bboxes))
        return batch


def _get_image_for_inference(sdata: SpatialData, image_key: str | None = None) -> DataArray | DataTree:
    if image_key is not None:
        return get_spatial_element(sdata.images, key=image_key)

    cell_image = sdata.attrs.get(SopaAttrs.CELL_SEGMENTATION)
    tissue_image = sdata.attrs.get(SopaAttrs.TISSUE_SEGMENTATION)

    assert cell_image is None or tissue_image is None or cell_image == tissue_image, (
        "When different images are existing for cell and tissue segmentation, you need to provide the `image_key` argument"
    )

    return get_spatial_element(sdata.images, key=cell_image or tissue_image)


def _get_extraction_parameters(
    image: DataArray | DataTree, level: int | None, magnification: int | None
) -> tuple[DataArray, int, float]:
    if isinstance(image, DataArray):
        assert level == 0, "Level must be 0 when using a DataArray"
        assert magnification is None, "Magnification must be None when the image is not multi scale"
        return image, 0, 1, 1

    if level < 0:
        assert isinstance(level, int) and level >= -len(image.keys()), "Invalid level"
        level = len(image.keys()) + level

    if magnification is None:
        tile_resize_factor, level_downsample = 1, 1
        if level is None:
            log.warning("Both level and magnification arguments are None. Using level=0 by default.")
            level = 0
    else:
        level, tile_resize_factor, level_downsample = _get_level_for_magnification(image, magnification)

    assert f"scale{level}" in image, f"Level {level} not found in the image scales"

    image = next(iter(image[f"scale{level}"].values()))

    return image, level, tile_resize_factor, level_downsample


def _get_level_for_magnification(image: DataArray | DataTree, magnification: int) -> int:
    """Return the best level for a given downsampling factor"""
    slide_metadata, backend = image.attrs.get("metadata", {}), image.attrs.get("backend")

    assert slide_metadata, "No `metadata` field found in the image attributes"
    assert backend is not None, "No backend found in the image metadata, can not infer the level for the magnification"
    assert "level_downsamples" in slide_metadata, "Missing `level_downsamples` in the image metadata"
    assert "properties" in slide_metadata, "Missing `properties` in the image metadata"

    objective_power = slide_metadata["properties"].get(f"{backend}.objective-power")
    mpp_x = slide_metadata["properties"].get(f"{backend}.mpp-x")

    if objective_power:
        downsample = int(objective_power) / magnification
    elif mpp_x:
        mpp_objective = min([80, 40, 20, 10, 5], key=lambda obj: abs(10 / obj - float(mpp_x)))
        downsample = mpp_objective / magnification
    else:
        raise ValueError("No objective-power or mpp-x information found in the metadata")

    if downsample < 1.0:
        log.warning(
            f"The requested magnification {magnification}x is higher than the objective power {objective_power}x. "
            f"Using the highest available magnification with upscaling."
        )

    level = _get_best_level_for_downsample(slide_metadata["level_downsamples"], downsample)
    level_downsample = slide_metadata["level_downsamples"][level]
    tile_resize_factor = level_downsample / downsample

    return level, tile_resize_factor, level_downsample


def _get_best_level_for_downsample(level_downsamples: list[float], downsample: float) -> int:
    """Return the best level for a given downsampling factor"""
    if downsample <= 1.0:
        return 0
    for level, ds in enumerate(level_downsamples):
        if ds > downsample:
            return level - 1
    return len(level_downsamples) - 1
