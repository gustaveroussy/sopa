import logging
from pathlib import Path

import xarray as xr
from spatialdata import SpatialData
from spatialdata.models import Image2DModel

from sopa.constants import SopaAttrs
from sopa.io.reader.utils import _default_image_kwargs, _image_int_dtype

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
    log.warning("The `aicsimageio` reader is deprecated. Use `sopa.io.bioio` instead.")
    image_models_kwargs, _ = _default_image_kwargs(image_models_kwargs, None)
    aics_kwargs = {} if aics_kwargs is None else aics_kwargs

    try:
        from aicsimageio import AICSImage
    except ImportError:
        raise ImportError("You need to install aicsimageio, e.g. by running `pip install aicsimageio`")

    xarr: xr.DataArray = AICSImage(path, **aics_kwargs).xarray_dask_data

    assert len(xarr.coords["T"]) == 1, f"Only one time dimension is supported, found {len(xarr.coords['T'])}."

    if len(xarr.coords["Z"]) > 1:
        log.info(f"3D image found, only reading {z_stack:=}")

    xarr = xarr.isel(T=0, Z=z_stack).rename({"C": "c", "Y": "y", "X": "x"})
    xarr = _image_int_dtype(xarr)

    image = Image2DModel.parse(xarr, c_coords=xarr.coords["c"].values, **image_models_kwargs)

    return SpatialData(images={"image": image}, attrs={SopaAttrs.CELL_SEGMENTATION: "image"})


def bioio(
    path: Path,
    z_stack: int = 0,
    image_models_kwargs: dict | None = None,
    bioio_kwargs: dict | None = None,
) -> SpatialData:
    """Read an image using [bioio](https://github.com/bioio-devs/bioio). It supports special formats such as `ND2`, `CZI`, `LIF`, or `DV`.

    !!! note "Extra dependencies"
        To use this reader, you'll need the `bioio` dependency (`pip install bioio`). You may need extra dependencies specific to your format, see their [documentation](https://bioio-devs.github.io/bioio/OVERVIEW.html#reader-installation).

    Args:
        path: Path to the image file
        z_stack: (Only for 3D images) Index of the stack in the z-axis to use.
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        bioio_kwargs: Keyword arguments passed to `bioio.BioImage`.

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs, _ = _default_image_kwargs(image_models_kwargs, None)
    bioio_kwargs = {} if bioio_kwargs is None else bioio_kwargs

    try:
        from bioio import BioImage
    except ImportError:
        raise ImportError("You need to install bioio, e.g. by running `pip install bioio`")

    xarr: xr.DataArray = BioImage(path, **bioio_kwargs).xarray_dask_data

    assert len(xarr.coords["T"]) == 1, f"Only one time dimension is supported, found {len(xarr.coords['T'])}."

    if len(xarr.coords["Z"]) > 1:
        log.info(f"3D image found, only reading {z_stack:=}")

    xarr = xarr.isel(T=0, Z=z_stack).rename({"C": "c", "Y": "y", "X": "x"})
    xarr = _image_int_dtype(xarr)

    image = Image2DModel.parse(xarr, c_coords=xarr.coords["c"].values, **image_models_kwargs)

    return SpatialData(images={"image": image}, attrs={SopaAttrs.CELL_SEGMENTATION: "image"})
