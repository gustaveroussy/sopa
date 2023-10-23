import logging
from pathlib import Path

from spatialdata import SpatialData

from .._constants import VALID_DIMENSIONS
from .._sdata import get_spatial_image

log = logging.getLogger(__name__)


def sanity_check(sdata: SpatialData, delete_table: bool = False):
    assert (
        len(sdata.images) > 0
    ), f"The spatialdata object has no image. Sopa is not designed for this."

    if len(sdata.images) > 1:
        log.warn(
            f"The spatialdata object has {len(sdata.images)} images. We advise to run sopa on one image (which can have multiple channels and multiple scales)"
        )

    if len(sdata.points) > 1:
        log.warn(
            f"The spatialdata object has {len(sdata.points)} points objects. It's easier to have only one (corresponding to transcripts), since sopa will use it directly without providing a key argument"
        )

    if len(sdata.points) == 0:
        log.warn("No transcripts found. Some tools from sopa will not be available.")

    image = get_spatial_image(sdata)
    assert (
        image.dims == VALID_DIMENSIONS
    ), f"Image must have the following three dimensions: {VALID_DIMENSIONS}. Found {image.dims}"

    if sdata.table is not None:
        if delete_table:
            log.info(
                "The table (i.e. `sdata.table`) will not be saved, since it will be created later by sopa"
            )
            del sdata.table
        else:
            raise ValueError(
                "sdata.table exists. Delete it you want to use sopa, to avoid conflicts with future table generation"
            )


def write_standardized(sdata: SpatialData, sdata_path: str, delete_table: bool = False):
    sanity_check(sdata, delete_table)

    log.info(f"Writing the following spatialdata object to {sdata_path}:\n{sdata}")

    sdata_path: Path = Path(sdata_path)
    if sdata_path.exists() and not any(sdata_path.iterdir()):
        sdata_path.rmdir()

    sdata.write(sdata_path)
