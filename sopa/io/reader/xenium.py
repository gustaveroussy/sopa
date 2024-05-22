# Updated from spatialdata-io: https://spatialdata.scverse.org/projects/io/en/latest/
# In the future, we will completely rely on spatialdata-io

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import dask.array as da
import packaging.version
from dask.dataframe import read_parquet
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel
from spatialdata.transformations import Identity, Scale
from spatialdata_io._constants._constants import XeniumKeys
from spatialdata_io.readers.xenium import _parse_version_of_xenium_analyzer

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
        - `experiment.xenium`: metadata file
        - `morphology_focus.ome.tif`: morphology image (or a directory, for recent versions of the Xenium)


    Args:
        path: Path to the Xenium directory containing all the experiment files
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object representing the Xenium experiment
    """
    path = Path(path)
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)
    image_models_kwargs["c_coords"] = [XeniumKeys.MORPHOLOGY_FOCUS_CHANNEL_0.value]

    with open(path / XeniumKeys.XENIUM_SPECS) as f:
        specs = json.load(f)

    points = {"transcripts": _get_points(path, specs)}

    with open(path / XeniumKeys.XENIUM_SPECS) as f:
        specs = json.load(f)
    # to trigger the warning if the version cannot be parsed
    version = _parse_version_of_xenium_analyzer(specs, hide_warning=False)

    images = {}
    if version is None or version < packaging.version.parse("2.0.0"):
        images["morphology_focus"] = _get_images(
            path,
            XeniumKeys.MORPHOLOGY_FOCUS_FILE,
            imread_kwargs,
            image_models_kwargs,
        )
    else:
        morphology_focus_dir = path / XeniumKeys.MORPHOLOGY_FOCUS_DIR
        files = [
            morphology_focus_dir / XeniumKeys.MORPHOLOGY_FOCUS_CHANNEL_IMAGE.value.format(i)
            for i in range(4)
        ]
        files = [f for f in files if f.exists()]
        if len(files) not in [1, 4]:
            raise ValueError(
                "Expected 1 (no segmentation kit) or 4 (segmentation kit) files in the morphology focus directory, "
                f"found {len(files)}: {files}"
            )

        if len(files) == 4:
            image_models_kwargs["c_coords"] = [
                XeniumKeys.MORPHOLOGY_FOCUS_CHANNEL_0,
                XeniumKeys.MORPHOLOGY_FOCUS_CHANNEL_1,
                XeniumKeys.MORPHOLOGY_FOCUS_CHANNEL_2,
                XeniumKeys.MORPHOLOGY_FOCUS_CHANNEL_3,
            ]

        images["morphology_focus"] = _get_images(
            morphology_focus_dir,
            files,
            imread_kwargs,
            image_models_kwargs,
        )

    return SpatialData(images=images, points=points)


def _get_points(path: Path, specs: dict[str, Any]):
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


def _get_images(
    path: Path,
    file: str | list[str],
    imread_kwargs: dict,
    image_models_kwargs: dict,
):
    if isinstance(file, list):
        image = da.concatenate([imread(f, **imread_kwargs) for f in file], axis=0)
    else:
        image = imread(path / file, **imread_kwargs)
    return Image2DModel.parse(
        image,
        transformations={"global": Identity()},
        dims=("c", "y", "x"),
        **image_models_kwargs,
    )
