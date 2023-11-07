import dask.array as da
import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Polygon, box

from sopa.segmentation import aggregate


def test_average_channels_geometries():
    image = np.random.randint(1, 10, size=(3, 8, 16))
    arr = da.from_array(image, chunks=(1, 8, 8))
    xarr = xr.DataArray(arr, dims=["c", "y", "x"])

    cell_size = 4
    cell_start = [(0, 0), (6, 2), (9, 3)]

    # One cell is on the first block, one is overlapping on both blocks, and one is on the last block
    cells = [box(x, y, x + cell_size - 1, y + cell_size - 1) for x, y in cell_start]

    means = aggregate._average_channels_geometries(xarr, cells)

    true_means = np.stack(
        [image[:, y : y + cell_size, x : x + cell_size].mean(axis=(1, 2)) for x, y in cell_start]
    )

    assert (means == true_means).all()


def test_get_cell_id():
    polygons = [box(10, 10, 20, 28), box(15, 18, 25, 22), box(30, 35, 34, 42)]
    gdf = gpd.GeoDataFrame(geometry=polygons)
    df = pd.DataFrame(
        {"x": [1.5, 16, 23, 67, 33, 19, 22, 10], "y": [15, 21, 34, 5, 40, 20, 21, 10]}
    )

    cell_id = aggregate._get_cell_id(gdf, df)

    assert list(cell_id) == [0, 1, 0, 0, 3, 1, 2, 1]


def test_count_transcripts():
    df_pandas = pd.DataFrame(
        {
            "x": [1, 2, 3, 7, 11, 1, 2, 2],
            "y": [1, 1, 2, 8, 0, 2, 3, 3],
            "gene": ["a", "a", "b", "c", "a", "c", "b", "b"],
        }
    )
    points = dd.from_pandas(df_pandas, npartitions=2)
    polygons = [
        Polygon(((1, 2), (3, 2), (3, 4), (1, 4))),
        Polygon(((0, 0), (2, 0), (2, 2), (0, 2))),
        Polygon(((0, 0), (3, 0), (3, 3), (0, 3))),
    ]

    gdf = gpd.GeoDataFrame(geometry=polygons)

    adata = aggregate._count_transcripts_geometries(gdf, points, "gene")
    expected = np.array([[0, 3, 1], [2, 0, 1], [2, 3, 1]])

    assert (adata.X.toarray() == expected).all()
