import numpy as np
import pytest
import spatialdata

import sopa


@pytest.mark.wsi_backends
def test_tiffslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="tiffslide")


@pytest.mark.wsi_backends
def test_openslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="openslide")


@pytest.mark.wsi_backends
def test_slideio_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="slideio")


@pytest.mark.wsi_backends
def test_region_matching():
    import numpy as np

    from sopa.io.reader._wsi_reader import get_reader

    wsi_tiffslide = get_reader("tiffslide")("tests/CMU-1-Small-Region.svs")
    wsi_openslide = get_reader("openslide")("tests/CMU-1-Small-Region.svs")
    wsi_slideio = get_reader("slideio")("tests/CMU-1-Small-Region.svs")
    wsi_xarray = get_reader("xarray")(sopa.io.wsi("tests/CMU-1-Small-Region.svs")["CMU-1-Small-Region"])

    location = (780, 660)
    size = (512, 556)

    region_tiffslide = wsi_tiffslide.read_region(location, level=0, size=size)
    region_openslide = wsi_openslide.read_region(location, level=0, size=size)
    region_slideio = wsi_slideio.read_region(location, level=0, size=size)
    region_xarray = wsi_xarray.read_region(location, level=0, size=size)

    assert (np.array(region_tiffslide) == np.array(region_openslide)).all(), (
        "Regions from tiffslide and openslide do not match"
    )
    assert (np.array(region_tiffslide) == np.array(region_slideio)).all(), (
        "Regions from tiffslide and slideio do not match"
    )
    assert (np.array(region_tiffslide) == np.array(region_xarray)).all(), (
        "Regions from tiffslide and xarray do not match"
    )


@pytest.mark.wsi_backends
@pytest.mark.parametrize("model", ["histo_ssl", "resnet50", "dinov2"])  # TODO: add conch and hoptimus0?
def test_deterministic_embedding(model: str):
    patch_width = 224

    sdata_on_disk = sopa.io.wsi("tests/CMU-1-Small-Region.svs")
    sdata_on_disk.write(f"tests/test_wsi_{model}.zarr")
    sdata_on_disk = spatialdata.read_zarr(f"tests/test_wsi_{model}.zarr")
    assert "backend" not in sdata_on_disk["CMU-1-Small-Region"].attrs  # fallback to xarray reader

    sdata_tiffslide = sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="tiffslide")

    sopa.patches.compute_embeddings(sdata_tiffslide, model, patch_width)

    for other_backend in ["slideio", "openslide"]:
        sdata_backend = sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend=other_backend)
        sopa.patches.compute_embeddings(sdata_backend, model, patch_width)

        assert np.isclose(sdata_tiffslide[f"{model}_embeddings"].X, sdata_backend[f"{model}_embeddings"].X).all()

    sopa.patches.compute_embeddings(sdata_on_disk, model, patch_width)

    assert np.isclose(sdata_tiffslide[f"{model}_embeddings"].X, sdata_on_disk[f"{model}_embeddings"].X).all()
