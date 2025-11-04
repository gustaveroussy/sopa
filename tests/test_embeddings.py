import pytest

import sopa


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
