import logging
import time
from typing import Callable

import dask as da
import numpy as np
import tqdm
from spatial_image import SpatialImage
from spatialdata import SpatialData, bounding_box_query
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale

from sopa.segmentation import Patches2D

log = logging.getLogger(__name__)


def _get_best_level_for_downsample(level_downsamples, downsample, epsilon=0):
    """return the best level for a given downsampling factor"""
    if downsample <= 1.0:
        return 0
    for lvl, ds in enumerate(level_downsamples):
        if ds > downsample + epsilon:
            return lvl - 1
    return level_count - 1


def _get_extraction_parameters(tiff_metadata, target_magnification, patch_width):
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
        log.error(f"Error retrieving the mpp for {{image_key}}, skipping tile extration.")
        return False

    lvl = _get_best_level_for_downsample(
        tiff_metadata["level_downsamples"], downsample, epsilon=0.01
    )
    resize_f = tiff_metadata["level_downsamples"][lvl] / downsample
    tile_s = int(patch_width * downsample)

    return lvl, resize_f, tile_s


def _get_patch(image: SpatialImage, box: list, lvl: int, resize_f: float):
    import cv2

    patch = bounding_box_query(image, ("y", "x"), box[:2][::-1], box[2:][::-1], "pixels")
    np_patch = np.array(patch[f"scale{lvl}"]["image"].transpose("y", "x", "c"))
    if resize_f != 1:
        shape = np_patch.shape
        dim = (int(shape[0] * resize_f), int(shape[1] * resize_f))
        np_patch = cv2.resize(np_patch, dim)
    return np_patch.transpose(2, 0, 1)


def embed_batch(model: str, device: str) -> Callable:

    import torch

    import sopa.embedding.models as models

    model = getattr(models, model)()
    model.eval().to(device)

    def _(patch: np.ndarray):
        torch_patch = torch.tensor(patch / 255.0, dtype=torch.float32)
        if len(torch_patch.shape) == 3:
            torch_patch = torch_patch.unsqueeze(0)
        with torch.no_grad():
            emb = model(torch_patch.to(device)).squeeze()
        return emb.cpu()

    return _


def patch_embedding(
    sdata: SpatialData,
    image_key: str | None,
    model_class: str,
    magnification: float | int,
    patch_width: float | int,
    patch_overlap: float | int = 0,
    batch_size: int = 32,
    num_workers: int = 1,
    device: str = "cpu",
):
    img = sdata.images[image_key]
    embedder = embed_batch(model=model_class, device=device)

    tiff_metadata = sdata.images[image_key].attrs["metadata"]
    lvl, resize_f, tile_s = _get_extraction_parameters(tiff_metadata, magnification, patch_width)

    slide_dim = tiff_metadata["level_dimensions"][0]
    # TODO: make this embedding size agnostic. At the moment it is not working for histoSSL
    output = np.zeros((1024, slide_dim[1] // tile_s, slide_dim[0] // tile_s), dtype="float32")

    patches = Patches2D(sdata, image_key, tile_s, 0)
    for i in tqdm.tqdm(range(0, len(patches), batch_size)):
        patch_boxes = patches[i : i + batch_size]
        get_batches = [da.delayed(_get_patch)(img, b, lvl, resize_f) for b in patch_boxes]
        batch = da.compute(*get_batches, num_workers=num_workers)
        emb = embedder(np.stack(batch))
        x = (np.array(patch_boxes) // tile_s)[:, 0]
        y = (np.array(patch_boxes) // tile_s)[:, 1]
        output[:, y, x] = emb.T

    feats_img = SpatialImage(output, dims=("c", "y", "x"))
    scale_image = Image2DModel.parse(
        feats_img, transformations={"pixels": Scale([tile_s, tile_s], axes=("x", "y"))}
    )
    scale_image.coords["y"] = tile_s * scale_image.coords["y"]
    scale_image.coords["x"] = tile_s * scale_image.coords["x"]

    sdata.add_image(model_class, scale_image)

    log.info(f"Tissue segmentation saved in sdata['model_class']")


if __name__ == "__main__":
    import spatialdata_plot

    from sopa.io import wsi
    from sopa.segmentation.tissue import hsv_otsu

    slide = wsi("CMU-1.svs")
    hsv_otsu(slide, "CMU-1")
    patch_embedding(slide, "CMU-1", "DINOv2Features", 10, 224, 0, 64, 4)

    slide.pl.render_images("DINOv2Features", channel=[1, 2, 3]).pl.show(save="lala.png")
