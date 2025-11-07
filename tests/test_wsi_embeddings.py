import numpy as np
import pytest
import spatialdata

import sopa
from sopa.patches._inference import Inference
from sopa.patches.infer import Patches2D, _get_image_for_inference


@pytest.mark.wsi
@pytest.mark.parametrize("model", ["resnet50", "dinov2"])
def test_deterministic_embedding(model: str):
    patch_width = 224

    sdata_on_disk = sopa.io.wsi("tests/CMU-1-Small-Region.svs")
    sdata_on_disk.write(f"tests/test_wsi_{model}.zarr")
    sdata_on_disk = spatialdata.read_zarr(f"tests/test_wsi_{model}.zarr")
    assert "backend" not in sdata_on_disk["CMU-1-Small-Region"].attrs  # fallback to xarray reader

    sdata_tiffslide = sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="tiffslide")

    sopa.patches.compute_embeddings(sdata_tiffslide, model, patch_width)
    sopa.patches.compute_embeddings(sdata_on_disk, model, patch_width)

    assert np.isclose(sdata_tiffslide[f"{model}_embeddings"].X, sdata_on_disk[f"{model}_embeddings"].X).all()


@pytest.mark.wsi
def test_resize_patches():
    sdata = sopa.io.wsi("tests/CMU-1-Small-Region.svs")
    image = _get_image_for_inference(sdata)

    patch_size = 253
    n_bboxes = 2

    infer = Inference(image, "dummy", patch_size, 0, 15, "cpu", False)
    patches = Patches2D(sdata, infer.image, infer.patch_width, 0)

    batch = infer._torch_batch(patches.bboxes[:n_bboxes])
    assert batch.shape == (n_bboxes, 3, patch_size, patch_size)
