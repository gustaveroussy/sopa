import pytest

import sopa


@pytest.mark.wsi
def test_openslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="openslide")


@pytest.mark.wsi
def test_region_matching():
    from sopa.patches.loader import get_reader

    image = sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="openslide", as_image=True)

    wsi_openslide = get_reader(image, "openslide")
    wsi_tiffslide = get_reader(image, "tiffslide")
    wsi_xarray = get_reader(image, "xarray")

    x, y = 780, 660
    width, height = 512, 556

    region_openslide = wsi_openslide.read_region(x, y, width, height, level=0)
    region_tiffslide = wsi_tiffslide.read_region(x, y, width, height, level=0)
    region_xarray = wsi_xarray.read_region(x, y, width, height, level=0)

    assert (region_openslide == region_xarray).all(), "Regions from openslide and xarray do not match"
    assert (region_tiffslide == region_xarray).all(), "Regions from tiffslide and xarray do not match"
