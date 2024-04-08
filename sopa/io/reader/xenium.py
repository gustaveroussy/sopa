# Updated from spatialdata-io: https://spatialdata.scverse.org/projects/io/en/latest/
# In the future, we will completely rely on spatialdata-io (when stable enough)

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from dask.dataframe import read_parquet
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel
from spatialdata.transformations import Identity, Scale
from spatialdata_io._constants._constants import XeniumKeys

from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def xenium(
    path: str | Path,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
) -> SpatialData:
    """Read Xenium data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.xenium.html).

    This function reads the following files:
        - `transcripts.parquet`: transcripts locations and names
        - `morphology_mip.ome.tif`: morphology image

    Args:
        path: Path to the Xenium directory containing all the experiment files
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object representing the Xenium experiment
    """
    path = Path(path)
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    with open(path / XeniumKeys.XENIUM_SPECS) as f:
        specs = json.load(f)

    points = {"transcripts": _get_points_xenium(path, specs)}

    images = {
        "morphology_mip": _get_images_xenium(
            path,
            XeniumKeys.MORPHOLOGY_MIP_FILE,
            imread_kwargs,
            image_models_kwargs,
        )
    }

    return SpatialData(images=images, points=points)


def _get_points_xenium(path: Path, specs: dict[str, Any]):
    table = read_parquet(path / XeniumKeys.TRANSCRIPTS_FILE)
    table["feature_name"] = table["feature_name"].apply(
        lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x),
        meta=("feature_name", "object"),
    )

    transform = Scale([1.0 / specs["pixel_size"], 1.0 / specs["pixel_size"]], axes=("x", "y"))
    points = PointsModel.parse(
        table,
        coordinates={
            "x": XeniumKeys.TRANSCRIPTS_X,
            "y": XeniumKeys.TRANSCRIPTS_Y,
            "z": XeniumKeys.TRANSCRIPTS_Z,
        },
        feature_key=XeniumKeys.FEATURE_NAME,
        instance_key=XeniumKeys.CELL_ID,
        transformations={"global": transform},
    )
    return points


def _get_images_xenium(
    path: Path,
    file: str,
    imread_kwargs: dict,
    image_models_kwargs: dict,
):
    image = imread(path / file, **imread_kwargs)
    return Image2DModel.parse(
        image,
        transformations={"global": Identity()},
        dims=("c", "y", "x"),
        c_coords=list(map(str, range(len(image)))),
        **image_models_kwargs,
    )
