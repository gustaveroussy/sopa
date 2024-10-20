from __future__ import annotations

import logging
from typing import Callable

import numpy as np
import torch
from datatree import DataTree
from spatialdata.transformations import Scale, Sequence, get_transformation
from xarray import DataArray

from . import models

log = logging.getLogger(__name__)


def _get_best_level_for_downsample(level_downsamples: list[float], downsample: float, epsilon: float = 0.01) -> int:
    """Return the best level for a given downsampling factor"""
    if downsample <= 1.0:
        return 0
    for level, ds in enumerate(level_downsamples):
        if ds > downsample + epsilon:
            return level - 1
    return len(level_downsamples) - 1


def _get_extraction_parameters(
    slide_metadata: dict,
    patch_width: int,
    level: int | None,
    magnification: int | None,
    backend: str,
) -> tuple[int, int, int, float, bool]:
    """
    Given the metadata for the slide, a target magnification and a patch width,
    it returns the best scale to get it from (level), a resize factor (resize_factor),
    the corresponding patch size at scale0 (patch_width) and downsample factor between
    scale0 and extraction level.
    """
    if level is None and magnification is None:
        log.warning("Both level and magnification arguments are None. Using level=0 by default.")
        level = 0

    if backend is None:
        log.warning("No backend found, using downsample=1")

    if magnification is None or backend is None:
        return level, 1, patch_width * 2**level, 1.0, True  # TODO: what if scaling != 2?

    if slide_metadata["properties"].get(f"{backend}.objective-power"):
        objective_power = int(slide_metadata["properties"].get(f"{backend}.objective-power"))
        downsample = objective_power / magnification

    elif slide_metadata["properties"].get(f"{backend}.mpp-x"):
        mppx = float(slide_metadata["properties"].get(f"{backend}.mpp-x"))

        mpp_objective = min([80, 40, 20, 10, 5], key=lambda obj: abs(10 / obj - mppx))
        downsample = mpp_objective / magnification
    else:
        return None, None, None, None, False

    level = _get_best_level_for_downsample(slide_metadata["level_downsamples"], downsample)
    resize_factor = slide_metadata["level_downsamples"][level] / downsample
    patch_width = int(patch_width * downsample)

    return level, resize_factor, patch_width, downsample, True


class Inference:
    def __init__(
        self,
        image: DataTree | DataArray,
        model: Callable | str,
        patch_width: int,
        level: int | None = 0,
        magnification: int | None = None,
        device: str = None,
    ):
        self.image = image
        self.patch_width = patch_width
        self.level = level
        self.magnification = magnification

        if isinstance(model, str):
            assert hasattr(
                models, model
            ), f"'{model}' is not a valid model name under `sopa.patches.models`. Valid names are: {', '.join(models.__all__)}"
            self.model_str = model
            self.model: torch.nn.Module = getattr(models, model)()
        else:
            self.model_str = model.__class__.__name__
            self.model: torch.nn.Module = model

        self.device = device
        if device:
            self.model.to(device)

        self._get_extraction_parameters()

        if isinstance(self.image, DataTree):
            self.image = next(iter(self.image[f"scale{self.level}"].values()))

    def _get_extraction_parameters(self):
        if isinstance(self.image, DataArray):
            self.resize_factor, self.patch_width_scale0, self.downsample = 1, self.patch_width, 1
            return

        slide_metadata = self.image.attrs.get("metadata", {})
        self.level, self.resize_factor, self.patch_width_scale0, self.downsample, success = _get_extraction_parameters(
            slide_metadata,
            self.patch_width,
            self.level,
            self.magnification,
            backend=self.image.attrs.get("backend"),
        )
        if not success:
            log.error("Error retrieving the image mpp, skipping tile embedding.")
            return False

    def _torch_resize(self, tensor: torch.Tensor, resize_factor: float):
        from torchvision.transforms import Resize

        dim = (int(tensor.shape[-2] * resize_factor), int(tensor.shape[-1] * resize_factor))
        return Resize(dim)(tensor)

    def _numpy_patch(
        self,
        box: tuple[int, int, int, int],
        patch_width: int,
    ) -> np.ndarray:
        """
        Extract a numpy patch from the image given a bounding box
        and pads a patch to a specific width since some patches might be smaller (e.g., on edges)
        """
        image_patch = self.image.sel(x=slice(box[0], box[2]), y=slice(box[1], box[3])).values

        pad_x, pad_y = patch_width - image_patch.shape[1], patch_width - image_patch.shape[2]
        return np.pad(image_patch, ((0, 0), (0, pad_x), (0, pad_y)))

    def _torch_batch(self, bboxes: np.ndarray):
        """Retrives a batch of patches using the bboxes"""
        extraction_patch_width = int(np.round(self.patch_width_scale0 / self.downsample / self.resize_factor))

        batch = np.array([self._numpy_patch(box, extraction_patch_width) for box in bboxes])
        batch = torch.tensor(batch, dtype=torch.float32) / 255.0

        if self.resize_factor != 1:
            batch = self._torch_resize(batch, self.resize_factor)

        return batch

    @torch.no_grad()
    def infer_bboxes(self, bboxes: np.ndarray) -> torch.Tensor:
        patches = self._torch_batch(bboxes)  # shape (B,3,Y,X)

        if len(patches.shape) == 3:
            patches = patches.unsqueeze(0)

        embedding = self.model(patches.to(self.device)).squeeze()
        return embedding.cpu()  # shape (B * output_dim)

    def get_patches_transformations(self, patch_overlap: float) -> dict[str, Sequence]:
        patch_step = self.patch_width_scale0 - self.downsample * patch_overlap
        to_image = Sequence([Scale([patch_step, patch_step], axes=("x", "y"))])
        image_transformations = get_transformation(self.image, get_all=True)
        return {cs: to_image.compose_with(t) for cs, t in image_transformations.items()}
