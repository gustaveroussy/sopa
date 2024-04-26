from __future__ import annotations

import logging
from pathlib import Path

import xarray as xr
from spatialdata import SpatialData
from spatialdata.models import Image2DModel

from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def aicsimageio(
    path: Path,
    z_stack: int = 0,
    image_models_kwargs: dict | None = None,
    aics_kwargs: dict | None = None,
) -> SpatialData:
    """Read an image using [AICSImageIO](https://github.com/AllenCellModeling/aicsimageio). It supports special formats such as `ND2`, `CZI`, `LIF`, or `DV`.

    !!! note "Extra dependencies"
        To use this reader, you'll need the `aicsimageio` dependency (`pip install aicsimageio`). To read `.czi` images, you'll also need to install `aicspylibczi` (for instance `pip install aicspylibczi`).

    Args:
        path: Path to the image file
        z_stack: (Only for 3D images) Index of the stack in the z-axis to use.
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        aics_kwargs: Keyword arguments passed to `aicsimageio.AICSImage`.

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs, _ = _default_image_kwargs(image_models_kwargs, None)
    aics_kwargs = {} if aics_kwargs is None else aics_kwargs

    try:
        from aicsimageio import AICSImage
    except ImportError:
        raise ImportError(
            "You need to install aicsimageio, e.g. by running `pip install aicsimageio`"
        )

    xarr: xr.DataArray = AICSImage(path, **aics_kwargs).xarray_dask_data

    assert (
        len(xarr.coords["T"]) == 1
    ), f"Only one time dimension is supported, found {len(xarr.coords['T'])}."

    if len(xarr.coords["Z"]) > 1:
        log.info(f"3D image found, only reading {z_stack:=}")

    xarr = xarr.isel(T=0, Z=z_stack).rename({"C": "c", "Y": "y", "X": "x"})

    image = Image2DModel.parse(xarr, c_coords=xarr.coords["c"].values, **image_models_kwargs)

    return SpatialData(images={"image": image})
