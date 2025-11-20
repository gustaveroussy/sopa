import numpy as np
import pytest
import spatialdata

import sopa


@pytest.mark.wsi
def test_patches_embeddings_inference_clustering():
    sdata = sopa.io.toy_dataset(length=200)

    with pytest.raises(AssertionError):  # two images, and no image_key provided
        sopa.patches.compute_embeddings(sdata, "dummy", 50)

    sopa.patches.compute_embeddings(sdata, "dummy", 50, image_key="he_image")

    assert sdata["dummy_embeddings"].shape == (4, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 25, level=-1, image_key="he_image")

    assert sdata["dummy_embeddings"].shape == (1, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 13, patch_overlap=3, level=-1, image_key="he_image")

    assert sdata["dummy_embeddings"].shape == (9, 3)

    sdata["he_image"].attrs["backend"] = "tiffslide"
    sdata["he_image"].attrs["metadata"] = {
        "level_downsamples": [1, 2, 4],
        "properties": {"tiffslide.objective-power": 50},
    }

    sopa.patches.compute_embeddings(sdata, "dummy", 10, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (9, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 11, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (4, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 50, magnification=100, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (16, 3)

    sopa.patches.cluster_embeddings(sdata, "dummy_embeddings")

    assert "cluster" in sdata["dummy_embeddings"].obs

    # Testing the xarray slicing for the read region
    sopa.settings.native_read_region = False
    sopa.patches.compute_embeddings(sdata, "dummy", 10, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (9, 3)


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
    from sopa.patches import TileLoader

    sdata = sopa.io.wsi("tests/CMU-1-Small-Region.svs")

    patch_size = 253
    n_bboxes = 2

    tile_loader = TileLoader(sdata, patch_size, magnification=15)

    batch = tile_loader[:n_bboxes]
    assert batch.shape == (n_bboxes, 3, patch_size, patch_size)
