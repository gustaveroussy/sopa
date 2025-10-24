import pytest

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
