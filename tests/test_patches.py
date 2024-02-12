import tempfile

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box
from spatialdata import SpatialData

from sopa._sdata import get_key
from sopa.segmentation.patching import Patches2D, _get_cell_id
from sopa.utils.data import uniform


@pytest.fixture
def sdata() -> SpatialData:
    sdata = uniform(length=512, cell_density=1e-3)
    return sdata


def test_patchify_image(sdata: SpatialData):
    image_key = get_key(sdata, "images")

    patches = Patches2D(sdata, image_key, 300, 100)
    assert len(patches) == 9

    patches = Patches2D(sdata, image_key, 512, 0)
    assert len(patches) == 1


def _patchify_transcripts(sdata: SpatialData, width: int, overlap: int) -> list[int]:
    with tempfile.TemporaryDirectory() as baysor_temp_dir:
        patches = Patches2D(sdata, "transcripts", width, overlap)
        return patches.patchify_transcripts(baysor_temp_dir)


def test_patchify_baysor(sdata: SpatialData):
    valid_indices = _patchify_transcripts(sdata, 30, 10)
    assert len(valid_indices) == 9

    valid_indices = _patchify_transcripts(sdata, 52, 0)
    assert len(valid_indices) == 1


def test_get_cell_id():
    polygons = [box(10, 10, 20, 28), box(15, 18, 25, 22), box(30, 35, 34, 42)]
    gdf = gpd.GeoDataFrame(geometry=polygons)
    df = pd.DataFrame(
        {"x": [1.5, 16, 23, 67, 33, 19, 22, 10], "y": [15, 21, 34, 5, 40, 20, 21, 10]}
    )

    cell_id = _get_cell_id(gdf, df)

    assert list(cell_id) == [0, 1, 0, 0, 3, 1, 2, 1]
