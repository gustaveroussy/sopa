from pathlib import Path
from typing import Literal

from spatialdata import SpatialData
from xarray import DataTree

from ...constants import SopaAttrs


def wsi(
    path: str | Path,
    as_image: bool = False,
    image_key: str = "wsi",
    backend: Literal["openslide", "tiffslide", "bioformats"] = "tiffslide",
) -> SpatialData | DataTree:
    """Read a WSI into a `SpatialData` object.

    !!! info "Internal reader"
        Internally, Sopa is using [wsidata.open_wsi](https://wsidata.readthedocs.io/en/latest/api/_autogen/wsidata.open_wsi.html) to open the image.

    Args:
        path: Path to the WSI
        as_image: If `True`, returns a, image instead of a `SpatialData` object
        image_key: The key (or name) to use for the image in the `SpatialData` object
        backend: The library to use as a backend in order to load the WSI. One of: `"openslide"`, `"tiffslide"`, `"bioformats"`. See the `reader` argument from [wsidata.open_wsi](https://wsidata.readthedocs.io/en/latest/api/_autogen/wsidata.open_wsi.html).

    Returns:
        A `SpatialData` object with a multiscale 2D-image of shape `(C, Y, X)`, or just the DataTree if `as_image=True`
    """
    try:
        from wsidata import open_wsi
    except ImportError:
        raise ImportError("Please install the wsi extra to use this function, e.g. via `pip install 'sopa[wsi]'`")

    sdata = open_wsi(
        path, reader=backend, attach_images=True, attach_thumbnail=False, image_key=image_key
    ).to_spatialdata()

    sdata.attrs[SopaAttrs.TISSUE_SEGMENTATION] = image_key

    metadata = sdata[image_key].attrs.copy()
    metadata["objective-power"] = _get_objective_power(metadata.get("raw", "{}"))

    sdata[image_key].attrs = {
        "metadata": metadata,
        "backend": backend,
        "path": str(path),
    }
    sdata[image_key].name = image_key

    if as_image:
        return sdata[image_key]

    return sdata


def _get_objective_power(raw_metadata) -> int | None:
    import ast

    try:
        return int(ast.literal_eval(raw_metadata)["objective-power"])
    except (KeyError, ValueError):
        return None
