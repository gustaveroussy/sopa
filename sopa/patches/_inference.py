import logging
from typing import Callable

import dask
import numpy as np
import torch
from xarray import DataArray, DataTree

from .. import settings
from ..io.reader._wsi_reader import get_reader
from . import models

log = logging.getLogger(__name__)


class Inference:
    def __init__(
        self,
        image: DataTree | DataArray,
        model: Callable | str,
        patch_width: int,
        level: int | None = 0,
        magnification: int | None = None,
        device: str | None = None,
        data_parallel: bool | list[int] = False,
    ):
        self.image, self.level, self.tile_resize_factor, self.level_downsample = _get_extraction_parameters(
            image, level, magnification
        )

        _backend = image.attrs.get("backend")
        _path = image.attrs.get("path")

        try:
            if settings.native_read_region and _backend is not None:
                self.slide = get_reader(_backend)(_path)
            else:
                self.slide = get_reader("xarray")(image)
        except Exception as e:
            log.warning(
                f"Exception raised for '{_backend}' and path '{_path}'. Falling back to xarray reader. Error: {e}"
            )
            self.slide = get_reader("xarray")(image)

        self.patch_width = int(patch_width / self.tile_resize_factor)
        self.resized_patch_width = patch_width

        self.model_str, self.model = self._instantiate_model(model)

        self.device = device
        if self.device is not None:
            self.model.to(device)

        if data_parallel:
            ids = data_parallel if isinstance(data_parallel, list) else list(range(torch.cuda.device_count()))
            self.model = torch.nn.DataParallel(self.model, device_ids=ids)

    def _instantiate_model(self, model: Callable | str) -> tuple[str, torch.nn.Module]:
        if isinstance(model, str):
            assert model in models.available_models, (
                f"'{model}' is not a valid model name. Valid names are: {', '.join(list(models.available_models.keys()))}"
            )

            return model, models.available_models[model]()

        return model.__class__.__name__, model

    def _torch_resize(self, tensor: torch.Tensor):
        from torchvision.transforms import Resize

        dim = (self.resized_patch_width, self.resized_patch_width)
        return Resize(dim)(tensor)

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
        return np.transpose(padded_patch, (2, 0, 1))

    def _torch_batch(self, bboxes: np.ndarray):
        """Retrives a batch of patches using the bboxes"""

        delayed_patches = [dask.delayed(self._numpy_patch)(box) for box in bboxes]
        batch = np.array(dask.compute(*delayed_patches))
        batch = torch.tensor(batch, dtype=torch.float32) / 255.0

        return batch if self.tile_resize_factor == 1 else self._torch_resize(batch)

    @torch.no_grad()
    def infer_bboxes(self, bboxes: np.ndarray) -> torch.Tensor:
        patches = self._torch_batch(bboxes)  # shape (B, C, Y, X)

        assert len(patches.shape) == 4
        embedding = self.model(patches.to(self.device))
        assert len(embedding.shape) == 2, "The model must have the following signature: (B, C, Y, X) -> (B, C)"

        return embedding.cpu()  # shape (B, output_dim)


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
