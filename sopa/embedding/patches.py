import logging
from typing import Callable

import dask as da
import numpy as np
import tqdm
from spatial_image import SpatialImage
from spatialdata import SpatialData, bounding_box_query
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale

from sopa._sdata import get_key
from sopa.segmentation import Patches2D

log = logging.getLogger(__name__)


def _get_best_level_for_downsample(
    level_downsamples: list, 
    downsample: float, 
    epsilon: float=0
):
    """return the best level for a given downsampling factor"""
    if downsample <= 1.0:
        return 0
    for level, ds in enumerate(level_downsamples):
        if ds > downsample + epsilon:
            return level - 1
    return len(level_downsamples) - 1


def _get_extraction_parameters(
    tiff_metadata: dict, 
    target_magnification: int, 
    patch_width: int
):
    """Given the metadata for the slide, a target magnification and a patch width
       it returns the best scale to get it from (lvl), a resize factor (resize_f) 
       and the corresponding patch size at scale0 (tile_s)"""
    if tiff_metadata["properties"].get("tiffslide.objective-power"):
        downsample = (
            int(tiff_metadata["properties"].get("tiffslide.objective-power")) / target_magnification
        )
    elif tiff_metadata["properties"].get("tiffslide.mpp-x"):
        obj2mpp = {80: 0.125, 40: 0.25, 20: 0.5, 10: 1.0, 5: 2.0}
        mppdiff = []
        for objpow, mpp in obj2mpp.items():
            mppdiff += [abs(mpp - float(tiff_metadata["properties"].get("tiffslide.mpp-x")))]
        idx = np.argmin([abs(mpp - 0.44177416504682804) for objpow, mpp in obj2mpp.items()])
        mpp_obj = list(obj2mpp.keys())[idx]
        downsample = mpp_obj / target_magnification
    else:
        log.error("Error retrieving the mpp for {{image_key}}, skipping tile extration.")
        return False

    level = _get_best_level_for_downsample(
        tiff_metadata["level_downsamples"], downsample, epsilon=0.01
    )
    resize_f = tiff_metadata["level_downsamples"][level] / downsample
    tile_s = int(patch_width * downsample)

    return level, resize_f, tile_s


def _get_patch(image: SpatialImage, box: list, level: int, resize_f: float):
    import cv2

    patch = bounding_box_query(image, ("y", "x"), box[:2][::-1], box[2:][::-1], "pixels")
    np_patch = np.array(patch[f"scale{level}"]["image"].transpose("y", "x", "c"))
    if resize_f != 1:
        shape = np_patch.shape
        dim = (int(shape[0] * resize_f), int(shape[1] * resize_f))
        np_patch = cv2.resize(np_patch, dim)
    return np_patch.transpose(2, 0, 1)


def embed_batch(model_name: str, device: str) -> Callable:
    import torch

    import sopa.embedding.models as models

    assert hasattr(
        models, model_name
    ), f"'{model_name}' is not a valid model name under `sopa.embedding.models`. Valid names are: {', '.join(models.__all__)}"

    model: torch.nn.Module = getattr(models, model_name)()
    model.eval().to(device)

    def _(patch: np.ndarray):
        torch_patch = torch.tensor(patch / 255.0, dtype=torch.float32)
        if len(torch_patch.shape) == 3:
            torch_patch = torch_patch.unsqueeze(0)
        with torch.no_grad():
            embedding = model(torch_patch.to(device)).squeeze()
        return embedding.cpu()

    return _


def patch_embedding(
    sdata: SpatialData,
    image_key: str | None,
    model_name: str,
    magnification: float | int,
    patch_width: float | int,
    patch_overlap: float | int = 0,
    batch_size: int = 32,
    num_workers: int = 1,
    device: str = "cpu",
):
    image_key = get_key(sdata, "images", image_key)
    image = sdata.images[image_key]
    tiff_metadata = image.attrs["metadata"]

    embedder = embed_batch(model_name=model_name, device=device)

    level, resize_f, tile_s = _get_extraction_parameters(tiff_metadata, magnification, patch_width)

    slide_dim = tiff_metadata["level_dimensions"][0]
    # TODO: make this embedding size agnostic. At the moment it is not working for histoSSL

    patches = Patches2D(sdata, image_key, tile_s, 0)
    output = np.zeros((1024, *patches.shape), dtype="float32")

    for i in tqdm.tqdm(range(0, len(patches), batch_size)):
        patch_boxes = patches[i : i + batch_size]

        get_batches = [da.delayed(_get_patch)(image, box, level, resize_f) for box in patch_boxes]
        batch = da.compute(*get_batches, num_workers=num_workers)
        embedding = embedder(np.stack(batch))

        xy = np.array([patches.pair_indices(k) for k in range(i, i + batch_size)]).T
        output[:, xy[1], xy[0]] = embedding.T

    embedding_image = SpatialImage(output, dims=("c", "y", "x"))
    embedding_image = Image2DModel.parse(
        embedding_image, transformations={"pixels": Scale([tile_s, tile_s], axes=("x", "y"))}
    )
    embedding_image.coords["y"] = tile_s * embedding_image.coords["y"]
    embedding_image.coords["x"] = tile_s * embedding_image.coords["x"]

    sdata.add_image(model_name, embedding_image)

    log.info(f"Tissue segmentation saved in sdata['{model_name}']")
