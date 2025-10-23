import sopa

def test_tiffslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="tiffslide")

def test_openslide_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="openslide")

def test_slideio_backend():
    sopa.io.wsi("tests/CMU-1-Small-Region.svs", backend="slideio")