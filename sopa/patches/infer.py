from __future__ import annotations

import logging
import warnings
from typing import Callable

import numpy as np
import tqdm
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from xarray import DataArray

from .._constants import SopaKeys
from ..segmentation import Patches2D
from ..utils import add_spatial_element, get_spatial_image

log = logging.getLogger(__name__)


def infer_wsi_patches(*args, **kwargs):
    warnings.warn(
        "`infer_wsi_patches` is deprecated, use `sopa.patches.compute_embeddings` instead",
        DeprecationWarning,
        stacklevel=2,
    )
    compute_embeddings(*args, **kwargs)


def compute_embeddings(
    sdata: SpatialData,
    model: Callable | str,
    patch_width: int,
    patch_overlap: int = 0,
    level: int | None = 0,
    magnification: int | None = None,
    image_key: str | None = None,
    batch_size: int = 32,
    device: str = None,
) -> DataArray:
    """Create an image made of patch based predictions of an image (useful for WSI images notably).

    !!! info
        The image will be saved into the `SpatialData` object with the key `sopa_{model_name}` (see the argument below).

    Args:
        sdata: A `SpatialData` object
        model: Callable that takes as an input a tensor of size (batch_size, channels, x, y) and returns a vector for each tile (batch_size, emb_dim), or a string with the name of one of the available models (`Resnet50Features`, `HistoSSLFeatures`, or `DINOv2Features`).
        patch_width: Width (pixels) of the patches.
        patch_overlap: Width (pixels) of the overlap between the patches.
        level: Image level on which the processing is performed. Either `level` or `magnification` should be provided.
        magnification: The target magnification on which the processing is performed. If `magnification` is provided, the `level` argument will be automatically computed.
        image_key: Optional image key of the image, unecessary if there is only one image.
        batch_size: Mini-batch size used during inference.
        device: Device used for the computer vision model.

    Returns:
        The `DataArray` of shape `(C,Y,X)` containing the model predictions (they are also saved in the `SpatialData` object).
    """
    try:
        import torch
    except ImportError:
        raise ImportError(
            "For patch embedding, you need `torch` (and perhaps `torchvision`). Consider installing the sopa WSI extra: `pip install 'sopa[wsi]'` (normal mode) or `pip install -e '.[wsi]'` (if using snakemake)"
        )

    from ._inference import Inference

    image_key, _ = get_spatial_image(sdata, key=image_key, return_key=True)
    image = sdata.images[image_key]

    infer = Inference(image, model, patch_width, level, magnification, device)
    patches = Patches2D(sdata, image_key, infer.patch_width_scale0, infer.downsample * patch_overlap)

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

    output_image = DataArray(output_image, dims=("c", "y", "x"))
    output_image = Image2DModel.parse(output_image, transformations=infer.get_patches_transformations(patch_overlap))

    output_key = f"sopa_{infer.model_str}"
    add_spatial_element(sdata, output_key, output_image)

    patches.write(shapes_key=SopaKeys.PATCHES_INFERENCE_KEY)

    return sdata[output_key]
