import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box
from spatialdata import SpatialData

import sopa
from sopa._constants import SopaKeys
from sopa.spatial.join import _get_cell_id
from sopa.utils.data import toy_dataset


@pytest.fixture
def sdata() -> SpatialData:
    sdata = toy_dataset(length=512, cell_density=1e-3)
    return sdata


def test_patchify_image(sdata: SpatialData):
    sopa.make_image_patches(sdata, 300, 100)
    assert len(sdata[SopaKeys.PATCHES]) == 9

    sopa.make_image_patches(sdata, 512, 0)
    assert len(sdata[SopaKeys.PATCHES]) == 1


def test_patchify_baysor(sdata: SpatialData):
    sopa.make_transcript_patches(sdata, 30, 10)
    assert len(sdata[SopaKeys.TRANSCRIPT_PATCHES]) == 9

    sopa.make_transcript_patches(sdata, 52, 0)
    assert len(sdata[SopaKeys.TRANSCRIPT_PATCHES]) == 1

    sopa.utils.delete_cache(sdata)


def test_get_cell_id():
    polygons = [box(10, 10, 20, 28), box(15, 18, 25, 22), box(30, 35, 34, 42)]
    gdf = gpd.GeoDataFrame(geometry=polygons)
    df = pd.DataFrame({"x": [1.5, 16, 23, 67, 33, 19, 22, 10], "y": [15, 21, 34, 5, 40, 20, 21, 10]})

    cell_id = _get_cell_id(gdf, df)

    assert list(cell_id) == [0, 1, 0, 0, 3, 1, 2, 1]
