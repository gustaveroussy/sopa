import dask_image.ndinterp
import numpy as np
from spatialdata import SpatialData
from xarray import DataArray, DataTree

from . import get_spatial_image


def resize_numpy(arr: np.ndarray, scale_factor: float, dims: list[str], output_shape: list[int]) -> np.ndarray:
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

    return dask_image.ndinterp.affine_transform(arr, matrix=transform, output_shape=output_shape).compute()


def assert_is_integer_dtype(dtype: np.dtype):
    assert np.issubdtype(dtype, np.integer), f"Expecting image to have an integer dtype, but found {dtype}"


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
    assert_is_integer_dtype(arr.dtype)
    assert_is_integer_dtype(dtype)

    if arr.dtype == dtype:
        return arr

    factor = np.iinfo(dtype).max / np.iinfo(arr.dtype).max
    return (arr * factor).astype(dtype)


def get_channel_names(image: DataArray | DataTree | SpatialData, image_key: str | None = None) -> np.ndarray:
    """Get the channel names of an image or a SpatialData object.

    Args:
        image: Either a `DataArray`, a `DataTree`, or a `SpatialData` object. If a `SpatialData` object, the `image_key` argument can be used.
        image_key: If `image` is a SpatialData object, the key of the image to get the channel names from. If `None`, tries to get it automatically.

    Returns:
        An array of channel names.
    """
    if isinstance(image, SpatialData):
        image = get_spatial_image(image, key=image_key)

    if isinstance(image, DataArray):
        return image.coords["c"].values
    if isinstance(image, DataTree):
        return image["scale0"].coords["c"].values
    raise ValueError(f"Image must be a DataTree or a DataArray. Found: {type(image)}")


def is_valid_c_coords(c_coords: np.ndarray) -> bool:
    return c_coords.dtype.kind in {"U", "S", "O"}


def ensure_string_channel_names(sdata: SpatialData, default_single_channel: str | None = "DAPI"):
    for key, image in list(sdata.images.items()):
        c_coords = get_channel_names(image)

        if is_valid_c_coords(c_coords):
            continue

        c_coords = [str(i) for i in range(len(c_coords))]
        if len(c_coords) == 1 and default_single_channel is not None:
            c_coords = [default_single_channel]

        sdata.set_channel_names(key, c_coords)
