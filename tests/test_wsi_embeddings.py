import numpy as np
import pytest
import spatialdata

import sopa


@pytest.mark.wsi
@pytest.mark.parametrize("model", ["histo_ssl", "dinov2"])
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
    from sopa.patches._inference import TileLoader
    from sopa.patches.infer import Patches2D, _get_image_for_inference

    sdata = sopa.io.wsi("tests/CMU-1-Small-Region.svs")
    image = _get_image_for_inference(sdata)

    patch_size = 253
    n_bboxes = 2

    tile_loader = TileLoader(image, patch_size, 0, 15)
    patches = Patches2D(sdata, tile_loader.image, tile_loader.patch_width, 0)

    batch = tile_loader.get_batch(patches.bboxes[:n_bboxes])
    assert batch.shape == (n_bboxes, 3, patch_size, patch_size)
