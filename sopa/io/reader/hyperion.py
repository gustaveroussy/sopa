from __future__ import annotations

import logging
from pathlib import Path

import dask.array as da
import numpy as np
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def hyperion(
    path: Path, image_models_kwargs: dict | None = None, imread_kwargs: dict | None = None
) -> SpatialData:
    """Read Hyperion data as a `SpatialData` object

    Args:
        path: Path to the directory containing the Hyperion `.tiff` images
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    files = [file for file in Path(path).iterdir() if file.suffix == ".tiff"]

    names = _get_channel_names_hyperion(files)
    image = da.concatenate(
        [imread(file, **imread_kwargs) for file in files],
        axis=0,
    )
    image = (image / image.max(axis=(1, 2)).compute()[:, None, None] * 255).astype(np.uint8)
    image = image.rechunk(chunks=image_models_kwargs["chunks"])

    log.info(f"Found channel names {names}")

    image_name = Path(path).absolute().stem
    image = Image2DModel.parse(
        image,
        dims=("c", "y", "x"),
        transformations={"pixels": Identity()},
        c_coords=names,
        **image_models_kwargs,
    )

    return SpatialData(images={image_name: image})


def _get_channel_names_hyperion(files: list[Path]):
    return [file.name[:-9].split("_")[1] for file in files]
