import dask.array as da
import dask_image.ndinterp
import numpy as np
import xarray as xr


def resize(xarr: xr.DataArray, scale_factor: float) -> da.Array:
    resize_dims = [dim in ["x", "y"] for dim in xarr.dims]
    transform = np.diag(
        [scale_factor if resize_dim else 1 for resize_dim in resize_dims]
    )
    output_shape = [
        size // scale_factor if resize_dim else size
        for size, resize_dim in zip(xarr.shape, resize_dims)
    ]

    return dask_image.ndinterp.affine_transform(
        xarr.data, matrix=transform, output_shape=output_shape
    )
