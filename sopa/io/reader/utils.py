from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Callable

import dask.array as da
import tifffile as tf
from dask_image.imread import imread
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

log = logging.getLogger(__name__)


def _default_image_kwargs(
    image_models_kwargs: dict | None = None, imread_kwargs: dict | None = None
) -> tuple[dict, dict]:
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs
    imread_kwargs = {} if imread_kwargs is None else imread_kwargs

    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 1024, 1024)

    if "scale_factors" not in image_models_kwargs:
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    return image_models_kwargs, imread_kwargs


def _deduplicate_names(df):
    is_duplicated = df[0].duplicated(keep=False)
    df.loc[is_duplicated, 0] += " (" + df.loc[is_duplicated, 1] + ")"
    return df[0].values


def _deduplicate_c_coords(c_coords: list[str]) -> list[str]:
    counter, res = defaultdict(int), []
    for channel in c_coords:
        if channel not in counter:
            res.append(channel)
        else:
            res.append(f"{channel} ({counter[channel]})")
        counter[channel] += 1
    return res


def _get_files_stem(files: list[Path]):
    return [file.stem for file in files]


def _general_tif_directory_reader(
    path: str,
    files_to_channels: Callable = _get_files_stem,
    suffix: str = ".tif",
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
):
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    files = [file for file in Path(path).iterdir() if file.suffix == suffix]

    names = files_to_channels(files)
    image = da.concatenate(
        [imread(file, **imread_kwargs) for file in files],
        axis=0,
    )
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


def _ome_channels_names(path: str):
    import xml.etree.ElementTree as ET

    tiff = tf.TiffFile(path)
    omexml_string = tiff.pages[0].description

    root = ET.fromstring(omexml_string)
    namespaces = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    channels = root.findall("ome:Image[1]/ome:Pixels/ome:Channel", namespaces)
    return [c.attrib["Name"] if "Name" in c.attrib else c.attrib["ID"] for c in channels]


def ome_tif(path: Path, as_image: bool = False) -> SpatialImage | SpatialData:
    """Read an `.ome.tif` image. This image should be a 2D image (with possibly multiple channels).
    Typically, this function can be used to open Xenium IF images.

    Args:
        path: Path to the `.ome.tif` image
        as_image: If `True`, will return a `SpatialImage` object

    Returns:
        A `SpatialImage` or a `SpatialData` object
    """
    image_models_kwargs, _ = _default_image_kwargs()
    image_name = Path(path).absolute().name.split(".")[0]
    image: da.Array = imread(path)

    if image.ndim == 4:
        assert image.shape[0] == 1, "4D images not supported"
        image = da.moveaxis(image[0], 2, 0)
        log.info(f"Transformed 4D image into a 3D image of shape (c, y, x) = {image.shape}")
    elif image.ndim != 3:
        raise ValueError(f"Number of dimensions not supported: {image.ndim}")

    image = image.rechunk(chunks=image_models_kwargs["chunks"])

    channel_names = _ome_channels_names(path)
    if len(channel_names) != len(image):
        channel_names = [str(i) for i in range(len(image))]
        log.warn(f"Channel names couldn't be read. Using {channel_names} instead.")

    image = SpatialImage(image, dims=["c", "y", "x"], name=image_name, coords={"c": channel_names})

    if as_image:
        return image

    image = Image2DModel.parse(
        image,
        dims=("c", "y", "x"),
        c_coords=channel_names,
        transformations={"pixels": Identity()},
        **image_models_kwargs,
    )

    return SpatialData(images={image_name: image})
