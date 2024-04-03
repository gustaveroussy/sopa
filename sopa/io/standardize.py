from __future__ import annotations

import logging
from pathlib import Path

import spatialdata
from spatialdata import SpatialData

from .._constants import VALID_DIMENSIONS, SopaKeys
from .._sdata import get_spatial_image
from ..utils.image import _check_integer_dtype

log = logging.getLogger(__name__)


def sanity_check(sdata: SpatialData, delete_table: bool = False, warn: bool = False):
    assert (
        len(sdata.images) > 0
    ), "The spatialdata object has no image. Sopa is not designed for this."

    if len(sdata.images) != 1:
        message = f"The spatialdata object has {len(sdata.images)} images. We advise to run sopa on one image (which can have multiple channels and multiple scales)"
        if warn:
            log.warn(message)
        else:
            raise ValueError(message)
    else:
        image = get_spatial_image(sdata)
        assert (
            image.dims == VALID_DIMENSIONS
        ), f"Image must have the following three dimensions: {VALID_DIMENSIONS}. Found {image.dims}"
        _check_integer_dtype(image.dtype)

    if len(sdata.points) > 1:
        log.warn(
            f"The spatialdata object has {len(sdata.points)} points objects. It's easier to have only one (corresponding to transcripts), since sopa will use it directly without providing a key argument"
        )

    # TODO: see https://github.com/scverse/spatialdata/issues/402
    # image_channels: np.ndarray = image.coords["c"].values
    # if image_channels.dtype.type is not np.str_:
    #     log.warn(f"Channel names are not strings. Converting {image_channels} to string values.")
    #     sdata[image_key].data = sdata[image_key].assign_coords(c=image_channels.astype(str))

    if SopaKeys.TABLE in sdata.tables:
        if delete_table:
            log.info(
                f"The table `sdata.tables['{SopaKeys.TABLE}']` will not be saved, since it will be created later by sopa"
            )
            del sdata.tables[SopaKeys.TABLE]


def read_zarr_standardized(path: str, warn: bool = False) -> SpatialData:
    sdata = spatialdata.read_zarr(path)
    sanity_check(sdata, warn=warn)
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

    if len(sdata.points) == 0:
        log.warn("No transcripts found. Some tools from sopa will not be available.")

    log.info(f"Writing the following spatialdata object to {sdata_path}:\n{sdata}")

    sdata_path: Path = Path(sdata_path)

    _check_can_write_zarr(sdata_path)
    sdata.write(sdata_path)
