from __future__ import annotations

import logging
from typing import Callable

try:
    import torch
except ImportError:
    raise ImportError(
        "For patch embedding, you need `torch` (and perhaps `torchvision`). Consider installing the sopa WSI extra: `pip install 'sopa[wsi]'` (normal mode) or `pip install -e '.[wsi]'` (if using snakemake)"
    )

import numpy as np
import tqdm
from multiscale_spatial_image import MultiscaleSpatialImage
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale

from .._constants import SopaKeys
from .._sdata import get_intrinsic_cs, get_key, save_image
from ..segmentation import Patches2D
from . import models

log = logging.getLogger(__name__)


def _get_best_level_for_downsample(
    level_downsamples: list[float], downsample: float, epsilon: float = 0.01
) -> int:
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
        log.warn("Both level and magnification arguments are None. Using level=0 by default.")
        level = 0

    if backend is None:
        log.warn("No backend found, using downsample=1")

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
        image: MultiscaleSpatialImage | SpatialImage,
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
            ), f"'{model}' is not a valid model name under `sopa.embedding.models`. Valid names are: {', '.join(models.__all__)}"
            self.model_str = model
            self.model: torch.nn.Module = getattr(models, model)()
        else:
            self.model_str = model.__class__.__name__
            self.model: torch.nn.Module = model

        self.device = device
        if device:
            self.model.to(device)

        self.cs = get_intrinsic_cs(None, image)

        self._get_extraction_parameters()

        if isinstance(self.image, MultiscaleSpatialImage):
            self.image = next(iter(self.image[f"scale{self.level}"].values()))

    def _get_extraction_parameters(self):
        if isinstance(self.image, SpatialImage):
            self.resize_factor, self.patch_width_scale0, self.downsample = 1, self.patch_width, 1
            return

        slide_metadata = self.image.attrs.get("metadata", {})
        self.level, self.resize_factor, self.patch_width_scale0, self.downsample, success = (
            _get_extraction_parameters(
                slide_metadata,
                self.patch_width,
                self.level,
                self.magnification,
                backend=self.image.attrs.get("backend"),
            )
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
        Extract a numpy patch from the MultiscaleSpatialImage given a bounding box
        and pads a patch to a specific width since some patches might be smaller (e.g., on edges)
        """
        image_patch = self.image.sel(x=slice(box[0], box[2]), y=slice(box[1], box[3])).values

        pad_x, pad_y = patch_width - image_patch.shape[1], patch_width - image_patch.shape[2]
        return np.pad(image_patch, ((0, 0), (0, pad_x), (0, pad_y)))

    def _torch_batch(self, bboxes: np.ndarray):
        """Retrives a batch of patches using the bboxes"""
        extraction_patch_width = int(
            np.round(self.patch_width_scale0 / self.downsample / self.resize_factor)
        )

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


def infer_wsi_patches(
    sdata: SpatialData,
    model: Callable | str,
    patch_width: int,
    patch_overlap: int = 0,
    level: int | None = 0,
    magnification: int | None = None,
    image_key: str | None = None,
    batch_size: int = 32,
    device: str = None,
) -> SpatialImage | bool:
    """Create an image made of patch based predictions of a WSI image.

    !!! info
        The image will be saved into the `SpatialData` object with the key `sopa_{model_name}` (see the argument below).

    Args:
        sdata: A `SpatialData` object
        model: Callable that takes as an input a tensor of size (batch_size, channels, x, y) and returns a vector for each tile (batch_size, emb_dim), or a string with the name of one of the available models (`Resnet50Features`, `HistoSSLFeatures`, or `DINOv2Features`).
        patch_width: Width (pixels) of the patches.
        patch_overlap: Width (pixels) of the overlap between the patches.
        level: Image level on which the processing is performed. Either `level` or `magnification` should be provided.
        magnification: The target magnification on which the processing is performed. If `magnification` is provided, the `level` argument will be automatically computed.
        image_key: Optional image key of the WSI image, unecessary if there is only one image.
        batch_size: Mini-batch size used during inference.
        device: Device used for the computer vision model.

    Returns:
        If the processing was successful, returns the `SpatialImage` of shape `(C,Y,X)` containing the model predictions, else `False`
    """
    image_key = get_key(sdata, "images", image_key)
    image = sdata.images[image_key]

    infer = Inference(image, model, patch_width, level, magnification, device)
    patches = Patches2D(
        sdata, image_key, infer.patch_width_scale0, infer.downsample * patch_overlap
    )

    log.info(f"Processing {len(patches)} patches extracted from level {infer.level}")

    predictions = []
    for i in tqdm.tqdm(range(0, len(patches), batch_size)):
        prediction = infer.infer_bboxes(patches.bboxes[i : i + batch_size])
        predictions.extend(prediction)
    predictions = torch.stack(predictions)

    if len(predictions.shape) == 1:
        predictions = torch.unsqueeze(predictions, 1)

    output_image = np.zeros((predictions.shape[1], *patches.shape), dtype=np.float32)
    for (loc_x, loc_y), pred in zip(patches.ilocs, predictions):
        output_image[:, loc_y, loc_x] = pred

    patch_step = infer.patch_width_scale0 - infer.downsample * patch_overlap
    output_image = SpatialImage(output_image, dims=("c", "y", "x"))
    output_image = Image2DModel.parse(
        output_image,
        transformations={infer.cs: Scale([patch_step, patch_step], axes=("x", "y"))},
    )

    output_key = f"sopa_{infer.model_str}"
    sdata.images[output_key] = output_image
    save_image(sdata, output_key)

    log.info(f"Patch predictions saved as an image in sdata['{output_key}']")

    patches.write(shapes_key=SopaKeys.PATCHES_INFERENCE_KEY)

    return sdata[output_key]
