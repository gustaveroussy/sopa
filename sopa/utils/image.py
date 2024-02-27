from __future__ import annotations

import dask.array as da
import dask_image.ndinterp
import numpy as np
import xarray as xr


def resize(xarr: xr.DataArray, scale_factor: float) -> da.Array:
    """Resize a xarray image

    Args:
        xarr: A `xarray` array
        scale_factor: Scale factor of resizing, e.g. `2` will decrease the width by 2

    Returns:
        Resized dask array
    """
    resize_dims = [dim in ["x", "y"] for dim in xarr.dims]
    transform = np.diag([scale_factor if resize_dim else 1 for resize_dim in resize_dims])
    output_shape = [
        size // scale_factor if resize_dim else size
        for size, resize_dim in zip(xarr.shape, resize_dims)
    ]

    return dask_image.ndinterp.affine_transform(
        xarr.data, matrix=transform, output_shape=output_shape
    )


def resize_numpy(
    arr: np.ndarray, scale_factor: float, dims: list[str], output_shape: list[int]
) -> np.ndarray:
    """Resize a numpy image

    Args:
        arr: a `numpy` array
        scale_factor: Scale factor of resizing, e.g. `2` will decrease the width by 2
        dims: List of dimension names. Only `"x"` and `"y"` are resized.
        output_shape: Size of the output array

    Returns:
        Resized array
    """
    resize_dims = [dim in ["x", "y"] for dim in dims]
    transform = np.diag([scale_factor if resize_dim else 1 for resize_dim in resize_dims])

    return dask_image.ndinterp.affine_transform(
        arr, matrix=transform, output_shape=output_shape
    ).compute()


def _check_integer_dtype(dtype: np.dtype):
    assert np.issubdtype(
        dtype, np.integer
    ), f"Expecting image to have an intenger dtype, but found {dtype}"


def scale_dtype(arr: np.ndarray, dtype: np.dtype) -> np.ndarray:
    """Change the dtype of an array but keep the scale compared to the type maximum value.

    !!! note "Example"
        For an array of dtype `uint8` being transformed to `np.uint16`, the value `255` will become `65535`

    Args:
        arr: A `numpy` array
        dtype: Target `numpy` data type

    Returns:
        A scaled `numpy` array with the dtype provided.
    """
    _check_integer_dtype(arr.dtype)
    _check_integer_dtype(dtype)

    if arr.dtype == dtype:
        return arr

    factor = np.iinfo(dtype).max / np.iinfo(arr.dtype).max
    return (arr * factor).astype(dtype)
