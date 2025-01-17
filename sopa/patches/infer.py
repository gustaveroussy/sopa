import logging
from typing import Callable

import numpy as np
import tqdm
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from xarray import DataArray

from .._constants import SopaAttrs, SopaKeys
from ..utils import add_spatial_element, get_spatial_element
from . import Patches2D

log = logging.getLogger(__name__)


def infer_wsi_patches(*args, **kwargs):
    log.warning("`infer_wsi_patches` is deprecated, use `sopa.patches.compute_embeddings` instead")
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
    device: str | None = None,
    key_added: str | None = None,
) -> None:
    """It creates patches, runs a computer vision model on each patch, and store the embeddings of each all patches as an image. This is mostly useful for WSI images.

    !!! info
        The image will be saved into the `SpatialData` object with the key `"{model_name}_embeddings"` (see the `model_name` argument below), except if `key_added` is provided.
        The shapes of the patches will be saved with the key `"embeddings_patches"`.

    Args:
        sdata: A `SpatialData` object
        model: Callable that takes as an input a tensor of size `(batch_size, channels, x, y)` and returns a vector for each tile `(batch_size, emb_dim)`, or a string with the name of one of the available models (`resnet50`, `histo_ssl`, or `dinov2`).
        patch_width: Width (pixels) of the patches.
        patch_overlap: Width (pixels) of the overlap between the patches.
        level: Image level on which the processing is performed. Either `level` or `magnification` should be provided.
        magnification: The target magnification on which the processing is performed. If `magnification` is provided, the `level` argument will be automatically computed.
        image_key: Optional image key of the image, unecessary if there is only one image.
        batch_size: Mini-batch size used during inference.
        device: Device used for the computer vision model.
        key_added: Optional name of the spatial element that will be added (storing the embeddings).
    """
    try:
        import torch
    except ImportError:
        raise ImportError(
            "For patch embedding, you need `torch` (and perhaps `torchvision`). Consider installing the sopa WSI extra: `pip install 'sopa[wsi]'`."
        )

    from ._inference import Inference

    image = get_spatial_element(
        sdata.images,
        key=image_key or sdata.attrs.get(SopaAttrs.CELL_SEGMENTATION) or sdata.attrs.get(SopaAttrs.TISSUE_SEGMENTATION),
    )

    infer = Inference(image, model, patch_width, level, magnification, device)
    patches = Patches2D(sdata, infer.image, infer.patch_width, patch_overlap)

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

    key_added = key_added or f"{infer.model_str}_embeddings"
    add_spatial_element(sdata, key_added, output_image)

    patches.add_shapes(key_added=SopaKeys.EMBEDDINGS_PATCHES)
