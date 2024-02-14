import logging
from typing import Callable

import dask as da
import numpy as np
import tqdm
from multiscale_spatial_image import MultiscaleSpatialImage
from spatial_image import SpatialImage
from spatialdata import SpatialData, bounding_box_query
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale

from sopa._sdata import get_intrinsic_cs, get_key
from sopa.segmentation import Patches2D

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
    tiff_metadata: dict, magnification: int, patch_width: int
) -> tuple[int, int, int, bool]:
    """
    Given the metadata for the slide, a target magnification and a patch width,
    it returns the best scale to get it from (level), a resize factor (resize_factor)
    and the corresponding patch size at scale0 (patch_width)
    """
    if tiff_metadata["properties"].get("tiffslide.objective-power"):
        downsample = (
            int(tiff_metadata["properties"].get("tiffslide.objective-power")) / magnification
        )
    elif tiff_metadata["properties"].get("tiffslide.mpp-x"):
        obj2mpp = {80: 0.125, 40: 0.25, 20: 0.5, 10: 1.0, 5: 2.0}
        mppx = float(tiff_metadata["properties"].get("tiffslide.mpp-x"))

        mppdiff = [(key, abs(mpp - mppx)) for key, mpp in obj2mpp.items()]
        mpp_obj = min(mppdiff, key=lambda x: x[1])[0]
        downsample = mpp_obj / magnification
    else:
        return None, None, None, False

    level = _get_best_level_for_downsample(tiff_metadata["level_downsamples"], downsample)
    resize_factor = tiff_metadata["level_downsamples"][level] / downsample
    patch_width = int(patch_width * downsample)

    return level, resize_factor, patch_width, True


def _numpy_patch(
    image: MultiscaleSpatialImage,
    box: list,
    level: int,
    resize_factor: float,
    coordinate_system: str,
) -> np.ndarray:
    import cv2

    multiscale_patch = bounding_box_query(
        image, ("y", "x"), box[:2][::-1], box[2:][::-1], coordinate_system
    )
    patch = np.array(
        next(iter(multiscale_patch[f"scale{level}"].values())).transpose("y", "x", "c")
    )

    if resize_factor != 1:
        dim = (int(patch.shape[0] * resize_factor), int(patch.shape[1] * resize_factor))
        patch = cv2.resize(patch, dim)

    return patch.transpose(2, 0, 1)


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

    return _, model.output_dim


def embed_wsi_patches(
    sdata: SpatialData,
    model_name: str,
    magnification: int,
    patch_width: int,
    image_key: str | None = None,
    batch_size: int = 32,
    num_workers: int = 1,
    device: str = "cpu",
):
    """Create an image made of patch embeddings of a WSI image.

    !!! info
        The image will be saved into the `SpatialData` object with the key `model_name` (see the argument below).

    Args:
        sdata: A `SpatialData` object
        model_name: Name of the computer vision model to be used. One of `Resnet50Features`, `HistoSSLFeatures`, or `DINOv2Features`.
        magnification: The target magnification.
        patch_width: Width of the patches for which the embeddings will be computed.
        image_key: Optional image key of the WSI image, unecessary if there is only one image.
        batch_size: Mini-batch size used during inference.
        num_workers: Number of workers used to extract patches.
        device: Device used for the computer vision model.
    """
    image_key = get_key(sdata, "images", image_key)
    image = sdata.images[image_key]

    assert isinstance(
        image, MultiscaleSpatialImage
    ), "Only `MultiscaleSpatialImage` images are supported"

    tiff_metadata = image.attrs["metadata"]
    coordinate_system = get_intrinsic_cs(sdata, image)

    embedder, output_dim = embed_batch(model_name=model_name, device=device)

    level, resize_factor, patch_width, success = _get_extraction_parameters(
        tiff_metadata, magnification, patch_width
    )
    if not success:
        log.error(f"Error retrieving the mpp for {image_key}, skipping tile embedding.")
        return False

    patches = Patches2D(sdata, image_key, patch_width, 0)
    output = np.zeros((output_dim, *patches.shape), dtype=np.float32)

    log.info(f"Computing {len(patches)} embeddings at level {level}")

    for i in tqdm.tqdm(range(0, len(patches), batch_size)):
        patch_boxes = patches[i : i + batch_size]

        get_batches = [
            da.delayed(_numpy_patch)(image, box, level, resize_factor, coordinate_system)
            for box in patch_boxes
        ]
        batch = da.compute(*get_batches, num_workers=num_workers)
        embedding = embedder(np.stack(batch))

        xy = np.array([patches.pair_indices(k) for k in range(i, i + batch_size)]).T
        output[:, xy[1], xy[0]] = embedding.T

    embedding_image = SpatialImage(output, dims=("c", "y", "x"))
    embedding_image = Image2DModel.parse(
        embedding_image,
        transformations={coordinate_system: Scale([patch_width, patch_width], axes=("x", "y"))},
    )
    embedding_image.coords["y"] = patch_width * embedding_image.coords["y"]
    embedding_image.coords["x"] = patch_width * embedding_image.coords["x"]

    sdata.add_image(model_name, embedding_image)

    log.info(f"Tissue segmentation saved in sdata['{model_name}']")

if __name__ == '__main__':
    import spatialdata_plot
    from sopa.io import wsi
    from sopa.segmentation.tissue import hsv_otsu

    slide = wsi('CMU-1.svs')
    hsv_otsu(slide, "CMU-1")
    embed_wsi_patches(slide, "Resnet50Features", 10, 256, image_key="CMU-1", batch_size=64, num_workers=8,device='cuda')

    slide.pl.render_images("Resnet50Features", channel=[1, 2, 3]).pl.show(save='lala.png')