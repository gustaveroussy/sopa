from __future__ import annotations

import logging

import dask
import numpy as np
from spatialdata import SpatialData
from xarray import DataArray, DataTree

from .. import settings
from ..constants import SopaAttrs, SopaKeys
from ..utils import get_spatial_element
from .patches import Patches2D

log = logging.getLogger(__name__)


class TileLoader:
    image: DataArray
    level: int
    level_downsample: int | None
    tile_resize_factor: float

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
        self.slide = get_reader(multiscale_image)

        self.init_extraction_parameters(multiscale_image, level, magnification)

        self.resized_patch_width = patch_width
        self.raw_patch_width = int(patch_width / self.tile_resize_factor)

        self.patches = Patches2D(sdata, self.image, self.raw_patch_width, patch_overlap=patch_overlap, roi_key=roi_key)

    def _numpy_patch(self, box: tuple[int, int, int, int]) -> np.ndarray:
        """
        Extract a numpy patch from the image given a bounding box
        and pads a patch to a specific width since some patches might be smaller (e.g., on edges)
        """
        x, y = box[0], box[1]
        width, height = box[2] - box[0], box[3] - box[1]

        if not isinstance(self.slide, RegionReader):
            x, y = x * self.level_downsample, y * self.level_downsample

        return self.slide.read_region(x, y, width, height, self.level)

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

    def init_extraction_parameters(
        self, multiscale_image: DataArray | DataTree, level: int | None, magnification: int | None
    ) -> tuple[DataArray, int, float]:
        if isinstance(multiscale_image, DataArray):
            assert level == 0, "Level must be 0 when using a DataArray"
            assert magnification is None, "Magnification must be None when the image is not multi scale"
            self.image = multiscale_image
            self.level, self.level_downsample, self.tile_resize_factor = 0, 1, 1
            return

        slide_metadata: dict[str, float | str] = multiscale_image.attrs.get("metadata", {})

        if magnification is not None:
            assert slide_metadata, "No `metadata` field found in the image attributes"
            assert "level_downsample" in slide_metadata, "Missing `level_downsample` in the image metadata"

            self.level, self.tile_resize_factor = _get_level_for_magnification(slide_metadata, magnification)
        else:
            self.tile_resize_factor = 1
            if level is None:
                log.warning("Both level and magnification arguments are None. Using level=0 by default.")
                self.level = 0
            else:
                if level < 0:
                    assert isinstance(level, int) and level >= -len(multiscale_image.keys()), "Invalid level"
                    level = len(multiscale_image.keys()) + level
                self.level = level

        assert f"scale{self.level}" in multiscale_image, f"Level {self.level} not found in the image scales"
        self.image = next(iter(multiscale_image[f"scale{self.level}"].values()))

        if not isinstance(self.slide, RegionReader):  # level_downsample is required
            assert "level_downsample" in slide_metadata, "Missing `level_downsample` in the image metadata"
            self.level_downsample = slide_metadata["level_downsample"][self.level]
        else:
            self.level_downsample = None


def _get_image_for_inference(sdata: SpatialData, image_key: str | None = None) -> DataArray | DataTree:
    if image_key is not None:
        return get_spatial_element(sdata.images, key=image_key)

    cell_image = sdata.attrs.get(SopaAttrs.CELL_SEGMENTATION)
    tissue_image = sdata.attrs.get(SopaAttrs.TISSUE_SEGMENTATION)

    assert cell_image is None or tissue_image is None or cell_image == tissue_image, (
        "When different images are existing for cell and tissue segmentation, you need to provide the `image_key` argument"
    )

    return get_spatial_element(sdata.images, key=cell_image or tissue_image)


def _get_level_for_magnification(slide_metadata: dict[str, float | str], magnification: int) -> int:
    """Return the best level for a given downsampling factor"""
    objective_power = slide_metadata.get("objective-power")
    mpp = slide_metadata.get("mpp")

    if objective_power:
        downsample = int(objective_power) / magnification
    elif mpp:
        mpp_objective = min([80, 40, 20, 10, 5], key=lambda obj: abs(10 / obj - float(mpp)))
        downsample = mpp_objective / magnification
    else:
        raise ValueError("No objective-power or mpp-x information found in the metadata")

    if downsample < 1.0:
        log.warning(
            f"The requested magnification {magnification}x is higher than the objective power {objective_power}x. "
            f"Using the highest available magnification with upscaling."
        )

    level = _get_best_level_for_downsample(slide_metadata["level_downsample"], downsample)
    tile_resize_factor = slide_metadata["level_downsample"][level] / downsample

    log.info(f"Using level {level} and tile_resize_factor {tile_resize_factor}")

    return level, tile_resize_factor


def _get_best_level_for_downsample(level_downsample: list[float], downsample: float) -> int:
    """Return the best level for a given downsampling factor"""
    if downsample <= 1.0:
        return 0
    for level, ds in enumerate(level_downsample):
        if ds > downsample:
            return level - 1
    return len(level_downsample) - 1


def get_reader(multiscale_image: DataTree | DataArray, name: str | None = None) -> RegionReader:
    name = name or (multiscale_image.attrs.get("backend") if settings.native_read_region else None)

    if name is None or name == "xarray":
        return RegionReader(multiscale_image)

    try:
        from wsidata import open_wsi
    except ImportError:
        raise ImportError("Please install the wsi extra to use this function, e.g. via `pip install 'sopa[wsi]'`")

    path = multiscale_image.attrs.get("path")
    assert path is not None, f"Found backend '{name}' but no path found in the image attributes"

    return open_wsi(path, reader=name)


class RegionReader:
    """Region-reader for xarray data (not using any specific backend like openslide or others)."""

    def __init__(self, multiscale_image: DataArray | DataTree):
        self.multiscale_image = multiscale_image

    def read_region(self, x: int, y: int, width: int, height: int, level: int = 0) -> np.ndarray:
        """Get a region from the slide."""
        image = (
            self.multiscale_image[f"scale{level}"].image
            if isinstance(self.multiscale_image, DataTree)
            else self.multiscale_image
        )
        tile = image[:, slice(y, y + height), slice(x, x + width)]
        patch: np.ndarray = tile.transpose("y", "x", "c").data.compute()

        pad_y, pad_x = height - patch.shape[0], width - patch.shape[1]
        padded_patch = np.pad(patch, ((0, pad_y), (0, pad_x), (0, 0)))

        return padded_patch
