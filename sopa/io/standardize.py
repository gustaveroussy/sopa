import logging

import spatialdata
from spatialdata import SpatialData

log = logging.getLogger(__name__)


def write_standardized(sdata: SpatialData, sdata_path: str):
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

    if sdata.table is not None:
        del sdata.table

    log.info(f"Writing the following spatialdata object to {sdata_path}:\n{sdata}")

    # sdata.write(sdata_path)
