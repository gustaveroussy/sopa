import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import spatialdata
from anndata import AnnData
from shapely import box
from spatialdata.models import ShapesModel

import sopa
from sopa.constants import SopaKeys
from sopa.spatial import cells_to_groups, mean_distance, niches_geometry_stats, spatial_neighbors, vectorize_niches
from sopa.spatial.distance import _random_distance_likelihood
from sopa.spatial.join import _get_cell_id

NICHE_KEY = "niche"
CT_KEY = "cell_type"


@pytest.fixture
def adata() -> AnnData:
    # c _ c - d _ d
    # | / - - - \ |
    # c - - - - - d
    # a - - - - - b
    # | \ - - - / |
    # a _ a - b _ b
    coords = np.array([
        [0, 0],
        [2, 0],
        [0, 2],
        [3, 0],
        [5, 0],
        [5, 2],
        [0, 3],
        [0, 5],
        [2, 5],
        [3, 5],
        [5, 5],
        [5, 3],
    ])
    niches = ["a", "a", "a", "b", "b", "b", "c", "c", "c", "d", "d", "d"]
    cell_types = ["c1", "c1", "c1", "c1", "c2", "c2", "c3", "c3", "c3", "c3", "c3", "c3"]

    adata = AnnData(obs=pd.DataFrame({NICHE_KEY: niches, CT_KEY: cell_types}))
    adata.obsm["spatial"] = coords

    spatial_neighbors(adata, radius=[0, 2.9])
    return adata


def test_cells_to_groups(adata: AnnData):
    df_distances = cells_to_groups(adata, NICHE_KEY, key_added_prefix=None)

    expected_0_0 = [0, 2, 2, 4]
    assert (df_distances.iloc[0] == np.array(expected_0_0)).all()

    expected_2_0 = [0, 1, 2, 3]
    assert (df_distances.iloc[1] == np.array(expected_2_0)).all()


def test_mean_distance(adata: AnnData):
    df_pairwise_distances = mean_distance(adata, NICHE_KEY)

    expected_a = [0, 5 / 3, 5 / 3, 10 / 3]
    assert (df_pairwise_distances.iloc[0] == expected_a).all()

    assert (
        (
            mean_distance(adata, CT_KEY, NICHE_KEY, ignore_zeros=True)
            >= mean_distance(adata, CT_KEY, NICHE_KEY, ignore_zeros=False)
        )
        .all()
        .all()
    )


def test_mean_distance_corrected(adata: AnnData):
    assert _random_distance_likelihood(0.5, ignore_zeros=True) > _random_distance_likelihood(0.5, ignore_zeros=False)

    diff = mean_distance(adata, group_key=CT_KEY, correction=False) - mean_distance(
        adata, group_key=CT_KEY, correction=True
    )

    assert (diff.values[np.diag_indices_from(diff)] == 0).all()

    diff = mean_distance(adata, group_key=CT_KEY, correction=False, ignore_zeros=True) - mean_distance(
        adata, group_key=CT_KEY, correction=True, ignore_zeros=True
    )

    assert (diff.std(0) < 1e-6).all()

    diff = mean_distance(adata, group_key=NICHE_KEY, target_group_key=CT_KEY, correction=False) - mean_distance(
        adata, group_key=NICHE_KEY, target_group_key=CT_KEY, correction=True
    )

    assert (diff.std(0) < 1e-6).all()
    assert diff["c3"].iloc[0] < diff["c1"].iloc[0] < diff["c2"].iloc[0]

    diff = mean_distance(adata, group_key=CT_KEY, target_group_key=NICHE_KEY, correction=False) - mean_distance(
        adata, group_key=CT_KEY, target_group_key=NICHE_KEY, correction=True
    )

    assert (diff.std() < 1e-6).all()  # same niche proportions


def test_vectorize_niches(adata: AnnData):
    gdf = vectorize_niches(adata, NICHE_KEY, buffer=0)

    assert len(gdf) == 4

    assert (gdf[SopaKeys.GEOMETRY_AREA] == 2).all()  # base * height / 2,

    assert (gdf[SopaKeys.GEOMETRY_LENGTH] == 4 + np.sqrt(8)).all()  # 2 + 2 + sqrt(2**2 + 2**2)


def test_niches_geometry_stats(adata: AnnData):
    df_geometries_stats = niches_geometry_stats(adata, NICHE_KEY, aggregation=["min", "mean"], buffer=0)

    expected_a = [0, 1, 1, np.sqrt(18)] * 2  # sqrt(3**2 + 3**2)
    assert (df_geometries_stats.iloc[0, 4:] == np.array(expected_a)).all()


def test_get_cell_id():
    polygons = [box(10, 10, 20, 28), box(15, 18, 25, 22), box(30, 35, 34, 42)]
    gdf = gpd.GeoDataFrame(geometry=polygons)
    df = pd.DataFrame({"x": [1.5, 16, 23, 67, 33, 19, 22, 10], "y": [15, 21, 34, 5, 40, 20, 21, 10]})

    cell_id = _get_cell_id(gdf, df)

    assert list(cell_id) == [0, 1, 0, 0, 3, 1, 2, 1]

    cell_id = _get_cell_id(gdf, df, None)  # should be [nan, 0, nan, nan, 2, 0, 1, 0]

    assert cell_id.isna().sum() == 3
    assert list(cell_id.iloc[[1, 4, 5, 6, 7]].values) == [0, 2, 0, 1, 0]


def test_sjoin():
    gdf = gpd.GeoDataFrame(
        geometry=[
            box(0, 0, 1, 1),
            box(1, 1, 2, 2),
            box(7, 7, 8, 8),
            box(4, 4, 6, 6),
        ]
    )
    gdf.geometry = gdf.geometry.buffer(0.1)
    gdf = ShapesModel.parse(gdf)

    gdf2 = gpd.GeoDataFrame(
        geometry=[
            box(0, 0, 1, 5),
            box(3, 3, 5, 5),
            box(10, 10, 11, 11),
        ]
    )
    gdf2.geometry = gdf2.geometry.buffer(0.1)
    gdf2 = ShapesModel.parse(gdf2)

    sdata = spatialdata.SpatialData(shapes={"left": gdf, "right": gdf2})

    index_right = sopa.spatial.sjoin(sdata, "left", "right")["index_right"]

    assert np.array_equal(index_right, [0, 0, np.nan, 1], equal_nan=True)

    index_right = sopa.spatial.sjoin(sdata, "right", "left")["index_right"]

    assert np.array_equal(index_right, [0, 1, 3, np.nan], equal_nan=True)
