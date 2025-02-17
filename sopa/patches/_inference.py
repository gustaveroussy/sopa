import logging
from typing import Callable

import numpy as np
import torch
from xarray import DataArray, DataTree

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
    ):
        self.image, self.level, self.resize_factor = _get_extraction_parameters(image, level, magnification)

        self.patch_width = int(patch_width / self.resize_factor)
        self.resized_patch_width = patch_width

        self.model_str, self.model = self._instantiate_model(model)

        self.device = device
        if self.device is not None:
            self.model.to(device)

    def _instantiate_model(self, model: Callable | str) -> tuple[str, torch.nn.Module]:
        if isinstance(model, str):
            assert (
                model in models.available_models
            ), f"'{model}' is not a valid model name. Valid names are: {', '.join(list(models.available_models.keys()))}"

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
        image_patch = self.image[:, slice(box[1], box[3]), slice(box[0], box[2])]

        pad_x, pad_y = self.patch_width - image_patch.shape[1], self.patch_width - image_patch.shape[2]
        return np.pad(image_patch, ((0, 0), (0, pad_x), (0, pad_y)))

    def _torch_batch(self, bboxes: np.ndarray):
        """Retrives a batch of patches using the bboxes"""

        batch = np.array([self._numpy_patch(box) for box in bboxes])
        batch = torch.tensor(batch, dtype=torch.float32) / 255.0

        return batch if self.resize_factor == 1 else self._torch_resize(batch)

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
        return image, 0, 1

    if level < 0:
        assert isinstance(level, int) and level >= -len(image.keys()), "Invalid level"
        level = len(image.keys()) + level

    if magnification is None:
        resize_factor = 1
        if level is None:
            log.warning("Both level and magnification arguments are None. Using level=0 by default.")
            level = 0
    else:
        level, resize_factor = _get_level_for_magnification(image, magnification)

    image = next(iter(image[f"scale{level}"].values()))

    return image, level, resize_factor


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

    if downsample <= 1.0:
        return 0, 1

    level = _get_best_level_for_downsample(slide_metadata["level_downsamples"], downsample)
    resize_factor = slide_metadata["level_downsamples"][level] / downsample

    return level, resize_factor


def _get_best_level_for_downsample(level_downsamples: list[float], downsample: float) -> int:
    """Return the best level for a given downsampling factor"""
    for level, ds in enumerate(level_downsamples):
        if ds > downsample:
            return level - 1
    return len(level_downsamples) - 1
