from pathlib import Path
from typing import Any, Literal

import xarray
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity, Scale
from xarray import DataArray, Dataset, DataTree

from ...constants import SopaAttrs


def wsi(
    path: str | Path,
    chunks: tuple[int, int, int] = (3, 512, 512),
    as_image: bool = False,
    backend: Literal["tiffslide", "openslide", "slideio"] = "tiffslide",
) -> SpatialData | DataTree:
    """Read a WSI into a `SpatialData` object.

    !!! info "Supported backends"
        Multiple backends are supported to read WSIs. To use `openslide`, you will need to install `openslide-python` and `openslide-bin`.
        If you want to use `slideio`, you will need to install `slideio`. The `tiffslide` backend is supported out-of-the-box.

    Args:
        path: Path to the WSI
        chunks: Tuple representing the chunksize for the dimensions `(C, Y, X)`.
        as_image: If `True`, returns a, image instead of a `SpatialData` object
        backend: The library to use as a backend in order to load the WSI. One of: `"openslide"`, `"tiffslide"`, `"slideio"`.

    Returns:
        A `SpatialData` object with a multiscale 2D-image of shape `(C, Y, X)`, or just the DataTree if `as_image=True`
    """
    image_name, img, slide_metadata = _open_wsi(path, backend=backend)

    images = {}
    for level, key in enumerate(sorted(img.keys(), key=int)):
        suffix = key if key != "0" else ""

        scale_image = DataArray(
            img[key].transpose("S", f"Y{suffix}", f"X{suffix}"),
            dims=("c", "y", "x"),
        ).chunk(chunks)

        scale_factor = slide_metadata["level_downsamples"][level]

        scale_image = Image2DModel.parse(
            scale_image[:3, :, :],
            transformations={"global": _get_scale_transformation(scale_factor)},
            c_coords=("r", "g", "b"),
        )
        scale_image.coords["y"] = scale_factor * scale_image.coords["y"]
        scale_image.coords["x"] = scale_factor * scale_image.coords["x"]

        images[f"scale{key}"] = Dataset({"image": scale_image})

    multiscale_image = DataTree.from_dict(images)
    sdata = SpatialData(images={image_name: multiscale_image}, attrs={SopaAttrs.TISSUE_SEGMENTATION: image_name})
    sdata[image_name].attrs["metadata"] = slide_metadata
    sdata[image_name].attrs["backend"] = backend
    sdata[image_name].attrs["path"] = str(path)
    sdata[image_name].name = image_name

    if as_image:
        return multiscale_image

    return sdata


def _get_scale_transformation(scale_factor: float):
    if scale_factor == 1:
        return Identity()
    return Scale([scale_factor, scale_factor], axes=("x", "y"))


def wsi_autoscale(
    path: str | Path,
    image_model_kwargs: dict | None = None,
    backend: Literal["tiffslide", "openslide", "slideio"] = "tiffslide",
) -> SpatialData:
    """Read a WSI into a `SpatialData` object.

    Scales are generated automatically by `spatialdata` instead of using
    the default multiscales.

    Args:
        path: Path to the WSI
        image_model_kwargs: Kwargs provided to the `Image2DModel`
        backend: The library to use as a backend in order to load the WSI. One of: `"openslide"`, `"tiffslide"`, `"slideio"`.

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_model_kwargs = _default_image_models_kwargs(image_model_kwargs)

    image_name, img, tiff_metadata = _open_wsi(path, backend=backend)

    img = img.rename_dims({"S": "c", "Y": "y", "X": "x"})

    multiscale_image = Image2DModel.parse(
        img["0"].transpose("c", "y", "x"),
        c_coords=("r", "g", "b"),
        **image_model_kwargs,
    )
    multiscale_image.attrs["metadata"] = tiff_metadata
    multiscale_image.attrs["backend"] = backend

    return SpatialData(images={image_name: multiscale_image}, attrs={SopaAttrs.TISSUE_SEGMENTATION: image_name})


def _default_image_models_kwargs(image_models_kwargs: dict | None) -> dict:
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs

    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (3, 4096, 4096)

    if "scale_factors" not in image_models_kwargs:
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    return image_models_kwargs


def _open_wsi(
    path: str | Path, backend: Literal["tiffslide", "openslide", "slideio"] = "openslide"
) -> tuple[str, xarray.Dataset, Any, dict]:
    from ._wsi_reader import get_reader

    image_name = Path(path).stem

    reader = get_reader(backend)(path)

    zarr_store = reader.get_zarr_store()
    metadata = reader.get_metadata()

    zarr_img = xarray.open_zarr(zarr_store, consolidated=False, mask_and_scale=False)

    return image_name, zarr_img, metadata
