import dask.array as da
import numpy as np
import xarray as xr
from shapely.geometry import box

from sopa.segmentation import shapes


def test_average_intensity():
    image = np.random.randint(1, 10, size=(3, 8, 16))
    arr = da.from_array(image, chunks=(1, 8, 8))
    xarr = xr.DataArray(arr, dims=["c", "y", "x"])

    cell_size = 4
    cell_start = [(0, 0), (6, 2), (9, 3)]

    # One cell is on the first block, one is overlapping on both blocks, and one is on the last block
    cells = [box(x, y, x + cell_size - 1, y + cell_size - 1) for x, y in cell_start]

    means = shapes.average(xarr, cells)

    true_means = np.stack(
        [image[:, y : y + cell_size, x : x + cell_size].mean(axis=(1, 2)) for x, y in cell_start]
    )

    assert (means == true_means).all()
