import pytest
import requests
import spatialdata
from spatialdata import SpatialData
from sopa.io import wsi
from sopa.patches import compute_embeddings

SVS_URL = "https://openslide.cs.cmu.edu/download/openslide-testdata/Aperio/CMU-1-Small-Region.svs"

@pytest.fixture(scope="session")
def svs_file(tmp_path_factory):
    """Download the SVS file into a temp directory."""
    tmp_dir = tmp_path_factory.mktemp("svs_data")
    svs_path = tmp_dir / "CMU-1-Small-Region.svs"

    if not svs_path.exists():
        response = requests.get(SVS_URL)
        response.raise_for_status()
        svs_path.write_bytes(response.content)

    return str(svs_path)

@pytest.fixture(scope="session")
def zarr_file(svs_file, tmp_path_factory):
    """Write the ZARR file based on the downloaded SVS file."""
    tmp_dir = tmp_path_factory.mktemp("zarr_data")
    zarr_path = tmp_dir / "CMU-1-Small-Region.zarr"

    slide = wsi(str(svs_file), backend="openslide")
    slide.write(str(tmp_dir / "CMU-1-Small-Region.zarr"), overwrite=True)

    return str(zarr_path)

@pytest.fixture(scope="session")
def sdata_backend(svs_file):
    return wsi(svs_file, backend="openslide")

@pytest.fixture(scope="session")
def sdata_zarr(zarr_file):
    return spatialdata.read_zarr(zarr_file)

def test_load_openslide(svs_file):
    slide = wsi(svs_file, backend="openslide")
    assert isinstance(slide, SpatialData)

def test_load_tiffslide(svs_file):
    slide = wsi(svs_file, backend="tiffslide")
    assert isinstance(slide, SpatialData)

def test_load_slideio(svs_file):
    slide = wsi(svs_file, backend="slideio")
    assert isinstance(slide, SpatialData)

def test_compute_embeddings_backend(sdata_backend):
    compute_embeddings(
        sdata_backend, 
        "resnet50", 
        256,
        patch_overlap=0, 
        magnification=20,
        image_key="CMU-1-Small-Region",
        batch_size=256, 
        device="cuda"
    )
    print(sdata_backend)
    assert sdata_backend["resnet50_embeddings"].X.shape[0] == 108

def test_compute_embeddings_zarr(sdata_zarr):
    compute_embeddings(
        sdata_zarr, 
        "resnet50", 
        256,
        patch_overlap=0, 
        magnification=20,
        image_key="CMU-1-Small-Region",
        batch_size=256, 
        device="cuda"
    )
    print(sdata_zarr)
    assert sdata_zarr["resnet50_embeddings"].X.shape[0] == 108
