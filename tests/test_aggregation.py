import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import box

from sopa.segmentation import aggregate


def test_average_channels():
    image = np.random.randint(1, 10, size=(3, 8, 16))
    arr = da.from_array(image, chunks=(1, 8, 8))
    xarr = xr.DataArray(arr, dims=["c", "y", "x"])

    cell_size = 4
    cell_start = [(0, 0), (6, 2), (9, 3)]

    # One cell is on the first block, one is overlapping on both blocks, and one is on the last block
    cells = [box(x, y, x + cell_size - 1, y + cell_size - 1) for x, y in cell_start]

    means = aggregate._average_channels(xarr, cells)

    true_means = np.stack(
        [image[:, y : y + cell_size, x : x + cell_size].mean(axis=(1, 2)) for x, y in cell_start]
    )

    assert (means == true_means).all()


def test_get_cell_id():
    polygons = [box(10, 10, 20, 28), box(15, 18, 25, 22), box(30, 35, 34, 42)]
    df = pd.DataFrame(
        {"x": [1.5, 16, 23, 67, 33, 19, 22, 10], "y": [15, 21, 34, 5, 40, 20, 21, 10]}
    )

    cell_id = aggregate._get_cell_id(polygons, df)

    assert list(cell_id) == [0, 1, 0, 0, 3, 1, 2, 0]
