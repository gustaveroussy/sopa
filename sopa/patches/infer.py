from __future__ import annotations

import logging

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
from spatialdata import SpatialData, bounding_box_query
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
) -> tuple[int, int, int, bool]:
    """
    Given the metadata for the slide, a target magnification and a patch width,
    it returns the best scale to get it from (level), a resize factor (resize_factor), 
    the corresponding patch size at scale0 (patch_width) and downsample factor between
    scale0 and extraction level.
    """
    if level is None and magnification is None:
        log.warn("Both level and magnification arguments are None. Using level=0 by default.")
        level = 0

    if magnification is None:
        return level, 1, patch_width * 2**level, True  # TODO: what if scaling != 2?

    if slide_metadata["properties"].get(f"{backend}.objective-power"):
        objective_power = int(slide_metadata["properties"].get(f"{backend}.objective-power"))
        downsample = objective_power / magnification

    elif slide_metadata["properties"].get(f"{backend}.mpp-x"):
        mppx = float(slide_metadata["properties"].get(f"{backend}.mpp-x"))

        mpp_objective = min([80, 40, 20, 10, 5], key=lambda obj: abs(10 / obj - mppx))
        downsample = mpp_objective / magnification
    else:
        return None, None, None, False

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
        device: str = "cpu",
    ):
        self.image = image

        if isinstance(model, str):
            assert hasattr(
                models, model
            ), f"'{model}' is not a valid model name under `sopa.embedding.models`. Valid names are: {', '.join(models.__all__)}"
            self.model_str = model
            self.model: torch.nn.Module = getattr(models, model)()
        else:
            self.model: torch.nn.Module = model

        self.device = device

        self.cs = get_intrinsic_cs(None, image)

        slide_metadata = image.attrs.get("metadata", {})
        self.level, self.resize_factor, self.patch_width, self.downsample, success = _get_extraction_parameters(
            slide_metadata, patch_width, level, magnification, backend=image.attrs['backend']
        )
        if not success:
            log.error("Error retrieving the image mpp, skipping tile embedding.")
            return False

        self.model.eval().to(device)

    def _resize(self, patch: np.ndarray):
        import cv2

        patch = patch.transpose(1, 2, 0)
        dim = (
            int(patch.shape[0] * self.resize_factor),
            int(patch.shape[1] * self.resize_factor),
        )
        patch = cv2.resize(patch, dim)
        return patch.transpose(2, 0, 1)

    def _torch_patch(
        self,
        box: tuple[int, int, int, int],
    ) -> np.ndarray:
        """Extract a numpy patch from the MultiscaleSpatialImage given a bounding box"""
        image_patch = bounding_box_query(
            self.image, ("y", "x"), box[:2][::-1], box[2:][::-1], self.cs
        )

        if isinstance(self.image, MultiscaleSpatialImage):
            image_patch = next(iter(image_patch[f"scale{self.level}"].values()))

        patch = image_patch.compute().data

        if self.resize_factor != 1:
            patch = self._resize(patch)

        return torch.tensor(patch / 255.0, dtype=torch.float32)

    def _torch_batch(self, bboxes: np.ndarray):
        batch = [self._torch_patch(box) for box in bboxes]

        max_y = max(img.shape[1] for img in batch)
        max_x = max(img.shape[2] for img in batch)

        def _pad(patch: torch.Tensor, max_y: int, max_x: int) -> torch.Tensor:
            pad_x, pad_y = max_x - patch.shape[2], max_y - patch.shape[1]
            return torch.nn.functional.pad(patch, (0, pad_x, 0, pad_y), value=0)

        return torch.stack([_pad(patch, max_y, max_x) for patch in batch])

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
    patch_overlap: int,
    level: int | None = 0,
    magnification: int | None = None,
    image_key: str | None = None,
    batch_size: int = 32,
    device: str = "cpu",
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
    patches = Patches2D(sdata, image_key, infer.patch_width,  infer.downsample*patch_overlap)

    log.info(f"Processing {len(patches)} patches extracted from level {infer.level}")

    predictions = []
    for i in tqdm.tqdm(range(0, len(patches), batch_size)):
        prediction = infer.infer_bboxes(patches.bboxes[i : i + batch_size])
        predictions.extend(prediction)
    predictions = torch.stack(predictions)

    output_image = np.zeros((predictions.shape[1], *patches.shape), dtype=np.float32)
    for (loc_x, loc_y), pred in zip(patches.ilocs, predictions):
        output_image[:, loc_y, loc_x] = pred

    patch_step = (infer.patch_width - infer.downsample * patch_overlap)
    output_image = SpatialImage(output_image, dims=("c", "y", "x"))
    output_image = Image2DModel.parse(
        output_image,
        transformations={
            infer.cs: Scale([patch_step, patch_step], axes=("x", "y"))
        },
    )

    output_key = f"sopa_{infer.model.__class__.__name__}"
    sdata.images[output_key] = output_image
    save_image(sdata, output_key)

    log.info(f"Patch predictions saved as an image in sdata['{output_key}']")

    patches.write(shapes_key=SopaKeys.PATCHES_INFERENCE_KEY)

    return sdata[output_key]
