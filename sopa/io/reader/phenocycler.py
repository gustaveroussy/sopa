# Readers for multiplex-imaging technologies
# In the future, we will completely rely on spatialdata-io (when all these functions exist)

from __future__ import annotations

import logging
import re
from pathlib import Path

import dask.array as da
import pandas as pd
import tifffile as tf
from dask.delayed import delayed
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

from .utils import _deduplicate_names, _default_image_kwargs

log = logging.getLogger(__name__)


def phenocycler(
    path: str | Path, channels_renaming: dict | None = None, image_models_kwargs: dict | None = None
) -> SpatialData:
    """Read Phenocycler data as a `SpatialData` object

    Args:
        path: Path to a `.qptiff` file, or a `.tif` file (if exported from QuPath)
        channels_renaming: A dictionnary whose keys correspond to channels and values to their corresponding new name. Not all channels need to be renamed.
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    image_models_kwargs, _ = _default_image_kwargs(image_models_kwargs)

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
