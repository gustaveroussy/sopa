from __future__ import annotations

from pathlib import Path
from typing import Any

import xarray
from multiscale_spatial_image import MultiscaleSpatialImage
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity, Scale


def wsi(
    path: str | Path, chunks: tuple[int, int, int] = (3, 256, 256), as_image: bool = False
) -> SpatialData:
    """Read a WSI into a `SpatialData` object

    Args:
        path: Path to the WSI
        chunks: Tuple representing the chunksize for the dimensions `(C, Y, X)`.
        as_image: If `True`, returns a, image instead of a `SpatialData` object

    Returns:
        A `SpatialData` object with a multiscale 2D-image of shape `(C, Y, X)`
    """
    image_name, img, tiff, tiff_metadata = _open_wsi(path)

    images = {}
    for level, key in enumerate(list(img.keys())):
        suffix = key if key != "0" else ""

        scale_image = SpatialImage(
            img[key].transpose("S", f"Y{suffix}", f"X{suffix}"),
            dims=("c", "y", "x"),
        ).chunk(chunks)

        scale_factor = tiff.level_downsamples[level]

        scale_image = Image2DModel.parse(
            scale_image,
            transformations={"pixels": _get_scale_transformation(scale_factor)},
            c_coords=("r", "g", "b"),
        )
        scale_image.coords["y"] = scale_factor * scale_image.coords["y"]
        scale_image.coords["x"] = scale_factor * scale_image.coords["x"]

        images[f"scale{key}"] = scale_image

    multiscale_image = MultiscaleSpatialImage.from_dict(images)
    multiscale_image.attrs["metadata"] = tiff_metadata
    multiscale_image.name = image_name

    if as_image:
        return multiscale_image

    return SpatialData(images={image_name: multiscale_image})


def _get_scale_transformation(scale_factor: float):
    if scale_factor == 1:
        return Identity()
    return Scale([scale_factor, scale_factor], axes=("x", "y"))


def wsi_autoscale(path: str | Path, image_model_kwargs: dict | None = None) -> SpatialData:
    """Read a WSI into a `SpatialData` object.

    Scales are generated automatically by `spatialdata` instead of using
    the default multiscales.

    Args:
        path: Path to the WSI
        image_model_kwargs: Kwargs provided to the `Image2DModel`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_model_kwargs = _default_image_models_kwargs(image_model_kwargs)

    image_name, img, _, tiff_metadata = _open_wsi(path)

    img = img.rename_dims({"S": "c", "Y": "y", "X": "x"})

    multiscale_image = Image2DModel.parse(
        img["0"].transpose("c", "y", "x"),
        transformations={"pixels": Identity()},
        c_coords=("r", "g", "b"),
        **image_model_kwargs,
    )
    multiscale_image.attrs["metadata"] = tiff_metadata

    return SpatialData(images={image_name: multiscale_image})


def _default_image_models_kwargs(image_models_kwargs: dict | None) -> dict:
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs

    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (3, 4096, 4096)

    if "scale_factors" not in image_models_kwargs:
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    return image_models_kwargs


def _open_wsi(path: str | Path) -> tuple[str, xarray.Dataset, Any, dict]:
    import tiffslide

    image_name = Path(path).stem

    tiff = tiffslide.open_slide(path)
    img = xarray.open_zarr(
        tiff.zarr_group.store,
        consolidated=False,
        mask_and_scale=False,
    )

    tiff_metadata = {
        "properties": tiff.properties,
        "dimensions": tiff.dimensions,
        "level_count": tiff.level_count,
        "level_dimensions": tiff.level_dimensions,
        "level_downsamples": tiff.level_downsamples,
    }
    return image_name, img, tiff, tiff_metadata
