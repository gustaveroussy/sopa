import numpy as np
import pytest

import sopa


@pytest.mark.wsi
def test_openslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="openslide")


@pytest.mark.wsi
def test_region_matching():
    from sopa.patches.loader import get_reader

    wsi_openslide = get_reader("openslide")("tests/CMU-1-Small-Region.svs")
    wsi_xarray = get_reader("xarray")(sopa.io.wsi("tests/CMU-1-Small-Region.svs")["wsi"])

    location = (780, 660)
    size = (512, 556)

    region_openslide = wsi_openslide.read_region(location, level=0, size=size)
    region_xarray = wsi_xarray.read_region(location, level=0, size=size)

    assert (np.array(region_openslide) == np.array(region_xarray)).all(), (
        "Regions from openslide and xarray do not match"
    )
