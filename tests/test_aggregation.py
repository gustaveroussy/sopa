import dask
import dask.array as da
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Polygon, box

from sopa.aggregation.channels import _aggregate_channels_aligned
from sopa.aggregation.transcripts import _count_transcripts_aligned

dask.config.set({"dataframe.query-planning": False})
import dask.dataframe as dd  # noqa


def test_aggregate_channels_aligned():
    image = np.random.randint(1, 10, size=(3, 8, 16))
    arr = da.from_array(image, chunks=(1, 8, 8))
    xarr = xr.DataArray(arr, dims=["c", "y", "x"])

    cell_size = 4
    cell_start = [(0, 0), (6, 2), (9, 3)]

    # One cell is on the first block, one is overlapping on both blocks, and one is on the last block
    cells = [box(x, y, x + cell_size - 1, y + cell_size - 1) for x, y in cell_start]

    mean_intensities = _aggregate_channels_aligned(xarr, cells, "average")
    min_intensities = _aggregate_channels_aligned(xarr, cells, "min")
    max_intensities = _aggregate_channels_aligned(xarr, cells, "max")

    true_mean_intensities = np.stack(
        [image[:, y : y + cell_size, x : x + cell_size].mean(axis=(1, 2)) for x, y in cell_start]
    )
    true_min_intensities = np.stack(
        [image[:, y : y + cell_size, x : x + cell_size].min(axis=(1, 2)) for x, y in cell_start]
    )
    true_max_intensities = np.stack(
        [image[:, y : y + cell_size, x : x + cell_size].max(axis=(1, 2)) for x, y in cell_start]
    )

    assert (mean_intensities == true_mean_intensities).all()
    assert (min_intensities == true_min_intensities).all()
    assert (max_intensities == true_max_intensities).all()


def test_count_transcripts():
    df_pandas = pd.DataFrame(
        {
            "x": [1, 2, 3, 7, 11, 1, 2, 2],
            "y": [1, 1, 2, 8, 0, 2, 3, 3],
            "gene": ["a", "a", "b", "c", "a", "c", "b", "b"],
        }
    )
    df_pandas["gene"] = df_pandas["gene"].astype(object)
    df_pandas["gene"].loc[0] = np.nan

    points = dd.from_pandas(df_pandas, npartitions=2)
    polygons = [
        Polygon(((1, 2), (3, 2), (3, 4), (1, 4))),
        Polygon(((0, 0), (2, 0), (2, 2), (0, 2))),
        Polygon(((0, 0), (3, 0), (3, 3), (0, 3))),
    ]

    gdf = gpd.GeoDataFrame(geometry=polygons)

    adata = _count_transcripts_aligned(gdf, points, "gene")
    expected = np.array([[0, 3, 1], [1, 0, 1], [1, 3, 1]])

    assert (adata.X.toarray() == expected).all()
