# Readers for multiplex-imaging technologies
# In the future, we will completely rely on spatialdata-io (when all these functions exist)


import logging
import re
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import tifffile as tf
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


def _default_image_models_kwargs(image_models_kwargs: dict | None):
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs

    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 4096, 4096)

    return image_models_kwargs


def macsima(
    path: Path, image_models_kwargs: dict | None = None, imread_kwargs: dict | None = None
) -> SpatialData:
    """Read MACSIMA data as a `SpatialData` object

    Notes:
        For all dulicated name, their index will be added in brackets after, for instance you will often find `DAPI (000)` to indicate the DAPI channel of index `000`

    Args:
        path: Path to the directory containing the MACSIMA `.tif` images
        image_models_kwargs: Kwargs provided to the `Image2DModel`
        imread_kwargs: Kwargs provided to `dask_image.imread.imread`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs = _default_image_models_kwargs(image_models_kwargs)
    imread_kwargs = {} if imread_kwargs is None else imread_kwargs

    files = [file for file in Path(path).iterdir() if file.suffix == ".tif"]

    names = _get_channel_names_macsima(files)
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


def qptiff(
    path: Path, channels_renaming: dict | None = None, image_models_kwargs: dict | None = None
) -> SpatialData:
    """Read a `.qptiff` file as a `SpatialData` object

    Args:
        path: Path to a `.qptiff` file
        channels_renaming: A dictionnary whose keys correspond to channels and values to their corresponding new name
        image_models_kwargs: Kwargs provided to the `Image2DModel`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs = _default_image_models_kwargs(image_models_kwargs)

    with tf.TiffFile(path) as tif:
        page_series = tif.series[0]
        names = _get_channel_names_qptiff(page_series)

        log.info(f"Found channel names {names}")

        if channels_renaming is not None and len(channels_renaming):
            log.info(f"Channels will be renamed by the dictionnary: {channels_renaming}")
            names = [channels_renaming.get(name, name) for name in names]
            log.info(f"New names are: {names}")

        image_name = Path(path).absolute().stem
        image = Image2DModel.parse(
            da.from_array(page_series.asarray(), chunks=image_models_kwargs["chunks"]),
            dims=list(page_series._axes.lower()),
            transformations={"pixels": Identity()},
            c_coords=names,
            **image_models_kwargs,
        )

        return SpatialData(images={image_name: image})


def phenocycler(path: Path, *args) -> SpatialData:
    """Read phenocycler data as a `SpatialData` object

    Args:
        path: Path to PhenoCycler `.qptiff` file
        args: Args provided to the `qptiff` reader function

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    return qptiff(path, *args)


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


def ome_tif(path: Path) -> SpatialImage:
    """Read an `.ome.tif` image. This image should be a 2D image (with possibly multiple channels).
    Typically, this function can be used to open Xenium IF images.

    Args:
        path: Path to the `.ome.tif` image

    Returns:
        A `SpatialImage`
    """
    image_models_kwargs = _default_image_models_kwargs(None)
    image_name = Path(path).absolute().name.split(".")[0]
    image: da.Array = imread(path)

    if image.ndim == 4:
        assert image.shape[0] == 1, f"4D images not supported"
        image = da.moveaxis(image[0], 2, 0)
        log.info(f"Transformed 4D image into a 3D image of shape (c, y, x) = {image.shape}")
    elif image.ndim != 3:
        raise ValueError(f"Number of dimensions not supported: {image.ndim}")

    image = image.rechunk(chunks=image_models_kwargs["chunks"])

    channel_names = _ome_channels_names(path)
    if len(channel_names) != len(image):
        channel_names = [str(i) for i in range(len(image))]
        log.warn(f"Channel names couldn't be read. Using {channel_names} instead.")

    return SpatialImage(image, dims=["c", "y", "x"], name=image_name, coords={"c": channel_names})
