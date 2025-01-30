import logging
from pathlib import Path

import spatialdata
from spatialdata import SpatialData

from .._constants import VALID_DIMENSIONS, SopaKeys
from ..utils import (
    assert_is_integer_dtype,
    get_channel_names,
    get_spatial_image,
    is_valid_c_coords,
)

log = logging.getLogger(__name__)


def sanity_check(sdata: SpatialData, delete_table: bool = False):
    assert len(sdata.images) > 0, "The spatialdata object has no image. Sopa is not designed for this."

    image = get_spatial_image(sdata)
    assert (
        image.dims == VALID_DIMENSIONS
    ), f"Image must have the following three dimensions: {VALID_DIMENSIONS}. Found {image.dims}"
    assert_is_integer_dtype(image.dtype)

    c_coords = get_channel_names(image)
    assert is_valid_c_coords(c_coords), f"Channel names must be strings, not {c_coords.dtype}"

    if SopaKeys.TABLE in sdata.tables:
        if delete_table:
            log.info(
                f"The table `sdata.tables['{SopaKeys.TABLE}']` will not be saved, since it will be created later by sopa"
            )
            del sdata.tables[SopaKeys.TABLE]


def read_zarr_standardized(path: str) -> SpatialData:
    sdata = spatialdata.read_zarr(path)
    sanity_check(sdata)
    return sdata


def _check_can_write_zarr(sdata_path: str):
    sdata_path = Path(sdata_path)

    if not sdata_path.exists():
        return

    assert not any(
        sdata_path.iterdir()
    ), f"Zarr directory {sdata_path} already exists. Sopa will not continue to avoid overwritting files."

    sdata_path.rmdir()


def write_standardized(sdata: SpatialData, sdata_path: str, delete_table: bool = False):
    sanity_check(sdata, delete_table)

    assert (
        SopaKeys.TABLE not in sdata.tables
    ), f"sdata.tables['{SopaKeys.TABLE}'] exists. Delete it you want to use sopa, to avoid conflicts with future table generation"

    log.info(f"Writing the following spatialdata object to {sdata_path}:\n{sdata}")

    sdata_path: Path = Path(sdata_path)

    _check_can_write_zarr(sdata_path)
    sdata.write(sdata_path)
