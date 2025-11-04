import numpy as np
import pytest

import sopa


@pytest.mark.wsi
def test_tiffslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="tiffslide")


@pytest.mark.wsi
def test_openslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="openslide")


@pytest.mark.wsi
def test_slideio_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="slideio")


@pytest.mark.wsi
def test_region_matching():
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
