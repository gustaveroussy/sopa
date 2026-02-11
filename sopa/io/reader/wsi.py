from pathlib import Path
from typing import Literal

from spatialdata import SpatialData
from xarray import DataTree

from ...constants import SopaAttrs


def wsi(
    path: str | Path,
    as_image: bool = False,
    backend: Literal["openslide", "bioformats"] = "openslide",
) -> SpatialData | DataTree:
    """Read a WSI into a `SpatialData` object.

    !!! info "Internal reader"
        Internally, Sopa is using [wsidata.open_wsi](https://wsidata.readthedocs.io/en/latest/api/_autogen/wsidata.open_wsi.html) to open the image.

    Args:
        path: Path to the WSI
        as_image: If `True`, returns a, image instead of a `SpatialData` object
        backend: The library to use as a backend in order to load the WSI. One of: `"openslide"`, `"bioformats"`. See the `reader` argument from [wsidata.open_wsi](https://wsidata.readthedocs.io/en/latest/api/_autogen/wsidata.open_wsi.html).

    Returns:
        A `SpatialData` object with a multiscale 2D-image of shape `(C, Y, X)`, or just the DataTree if `as_image=True`
    """
    try:
        from wsidata import open_wsi
    except ImportError:
        raise ImportError("Please install the wsi extra to use this function, e.g. via `pip install 'sopa[wsi]'`")

    sdata = open_wsi(path, reader=backend, attach_images=True, attach_thumbnail=False).to_spatialdata()

    image_name = "wsi"

    sdata.attrs[SopaAttrs.TISSUE_SEGMENTATION] = image_name

    metadata = sdata[image_name].attrs.copy()
    sdata[image_name].attrs = {
        "metadata": metadata,
        "backend": backend,
        "path": str(path),
    }
    sdata[image_name].name = image_name

    if as_image:
        return sdata[image_name]

    return sdata
