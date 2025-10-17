import logging

import spatialdata
from spatialdata import SpatialData

from .._constants import VALID_DIMENSIONS
from ..utils import (
    assert_is_integer_dtype,
    get_channel_names,
    get_spatial_image,
    is_valid_c_coords,
)

log = logging.getLogger(__name__)


def sanity_check(sdata: SpatialData):
    assert len(sdata.images) > 0, "The spatialdata object has no image. Sopa is not designed for this."

    image = get_spatial_image(sdata)
    assert image.dims == VALID_DIMENSIONS, (
        f"Image must have the following three dimensions: {VALID_DIMENSIONS}. Found {image.dims}"
    )
    assert_is_integer_dtype(image.dtype)

    c_coords = get_channel_names(image)
    assert is_valid_c_coords(c_coords), f"Channel names must be strings, not {c_coords.dtype}"


def read_zarr_standardized(path: str) -> SpatialData:
    sdata = spatialdata.read_zarr(path)

    sanity_check(sdata)

    return sdata
