# Readers for multiplex-imaging technologies
# In the future, we will completely rely on spatialdata-io (when all these functions exist)

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Callable

import dask.array as da
import numpy as np
import pandas as pd
import tifffile as tf
from dask.delayed import delayed
from dask_image.imread import imread
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

log = logging.getLogger(__name__)


def _deduplicate_names(df):
    is_duplicated = df[0].duplicated(keep=False)
    df.loc[is_duplicated, 0] += " (" + df.loc[is_duplicated, 1] + ")"
    return df[0].values


def _parse_name_macsima(file):
    index = file.name[2:5] if file.name[0] == "C" else file.name[:3]
    match = re.search(r"_A-(.*?)_C-", file.name)
    if match:
        antibody = match.group(1)
        channel = re.search(r"_C-(.*?)\.tif", file.name).group(1)
        uid = f"{channel}-{index}"
    else:
        antibody = re.search(r"_A-(.*?)\.tif", file.name).group(1)
        uid = index
    return [antibody, uid]


def _get_channel_names_macsima(files):
    df_antibodies = pd.DataFrame([_parse_name_macsima(file) for file in files])
    return _deduplicate_names(df_antibodies)


def _default_image_models_kwargs(image_models_kwargs: dict | None = None):
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs

    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 1024, 1024)

    if "scale_factors" not in image_models_kwargs:
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    return image_models_kwargs


def macsima(path: Path, **kwargs: int) -> SpatialData:
    """Read MACSIMA data as a `SpatialData` object

    Notes:
        For all dulicated name, their index will be added in brackets after, for instance you will often find `DAPI (000)` to indicate the DAPI channel of index `000`

    Args:
        path: Path to the directory containing the MACSIMA `.tif` images
        kwargs: Kwargs for `_general_tif_directory_reader`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    return _general_tif_directory_reader(
        path, files_to_channels=_get_channel_names_macsima, **kwargs
    )


def _get_files_stem(files: list[Path]):
    return [file.stem for file in files]


def _general_tif_directory_reader(
    path: str,
    files_to_channels: Callable = _get_files_stem,
    suffix: str = ".tif",
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
):
    image_models_kwargs = _default_image_models_kwargs(image_models_kwargs)
    imread_kwargs = {} if imread_kwargs is None else imread_kwargs

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


def _get_channel_name_qptiff(description):
    import xml.etree.ElementTree as ET

    root = ET.fromstring(description)

    for xml_path in [".//Biomarker", ".//ExcitationFilter//Bands//Name"]:
        field = root.find(xml_path)
        if field is not None:
            return field.text

    return re.search(r"<Name>(.*?)</Name>", description).group(1)


def _get_channel_names_qptiff(page_series):
    df_names = pd.DataFrame(
        [[_get_channel_name_qptiff(page.description), str(i)] for i, page in enumerate(page_series)]
    )
    return _deduplicate_names(df_names)


def _get_IJ_channel_names(path: str) -> list[str]:
    with tf.TiffFile(path) as tif:
        default_names = [str(i) for i in range(len(tif.pages))]

        if len(tif.pages) > 1:
            ij_metadata_tag = tif.pages[0].tags.get("IJMetadata", None)

            if ij_metadata_tag and "Labels" in ij_metadata_tag.value:
                return ij_metadata_tag.value["Labels"]

            log.warn("Could not find channel names in IJMetadata.")
            return default_names

        log.warn("The TIF file does not have multiple channels.")
        return default_names


def _rename_channels(names: list[str], channels_renaming: dict | None = None):
    log.info(f"Found channel names {names}")
    if channels_renaming is not None and len(channels_renaming):
        log.info(f"Channels will be renamed by the dictionnary: {channels_renaming}")
        names = [channels_renaming.get(name, name) for name in names]
        log.info(f"New names are: {names}")
    return names


def phenocycler(
    path: str | Path, channels_renaming: dict | None = None, image_models_kwargs: dict | None = None
) -> SpatialData:
    """Read Phenocycler data as a `SpatialData` object

    Args:
        path: Path to a `.qptiff` file, or a `.tif` file (if exported from QuPath)
        channels_renaming: A dictionnary whose keys correspond to channels and values to their corresponding new name. Not all channels need to be renamed.
        image_models_kwargs: Kwargs provided to the `Image2DModel`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs = _default_image_models_kwargs(image_models_kwargs)

    path = Path(path)
    image_name = path.absolute().stem

    if path.suffix == ".qptiff":
        with tf.TiffFile(path) as tif:
            series = tif.series[0]
            names = _get_channel_names_qptiff(series)

            delayed_image = delayed(lambda series: series.asarray())(tif)
            image = da.from_delayed(delayed_image, dtype=series.dtype, shape=series.shape)
    elif path.suffix == ".tif":
        image = imread(path)
        names = _get_IJ_channel_names(path)
    else:
        raise ValueError(f"Unsupported file extension {path.suffix}. Must be '.qptiff' or '.tif'.")

    names = _rename_channels(names, channels_renaming)
    image = image.rechunk(chunks=image_models_kwargs["chunks"])

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


def hyperion(
    path: Path, image_models_kwargs: dict | None = None, imread_kwargs: dict | None = None
) -> SpatialData:
    """Read Hyperion data as a `SpatialData` object

    Args:
        path: Path to the directory containing the Hyperion `.tiff` images
        image_models_kwargs: Kwargs provided to the `Image2DModel`
        imread_kwargs: Kwargs provided to `dask_image.imread.imread`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs = _default_image_models_kwargs(image_models_kwargs)
    imread_kwargs = {} if imread_kwargs is None else imread_kwargs

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
    image_models_kwargs = _default_image_models_kwargs()
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
