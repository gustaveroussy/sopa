import logging
import re
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import tifffile as tf
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

log = logging.getLogger(__name__)


def _deduplicate_names(df):
    is_duplicated = df[0].duplicated(keep=False)
    df.loc[is_duplicated, 0] += " (" + df.loc[is_duplicated, 1] + ")"
    return df[0].values


def _parse_name_macsima(file):
    index = file.name[:3]
    match = re.search(r"_A-(.*?)_C-", file.name)
    if match:
        antibody = match.group(1)
        channel = re.search(r"_C-(.*?)\.ome\.tif", file.name).group(1)
        uid = f"{channel}-{index}"
    else:
        antibody = re.search(r"_A-(.*?)\.ome\.tif", file.name).group(1)
        uid = index
    return [antibody, uid]


def _get_channel_names_macsima(files):
    df_antibodies = pd.DataFrame([_parse_name_macsima(file) for file in files])
    return _deduplicate_names(df_antibodies)


def macsima(
    path: Path, image_models_kwargs: dict | None = None, imread_kwargs: dict | None = None
) -> SpatialData:
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs
    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 4096, 4096)
    imread_kwargs = {} if imread_kwargs is None else imread_kwargs

    files = [file for file in Path(path).iterdir() if file.suffix == ".tif"]

    names = _get_channel_names_macsima(files)
    image = da.concatenate(
        [imread(file, **imread_kwargs) for file in files],
        axis=0,
    )

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
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs
    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 4096, 4096)

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


def phenocycler(path: Path, *args):
    return qptiff(path, *args)


def _get_channel_names_hyperion(files: list[Path]):
    return [file.name[:-9].split("_")[1] for file in files]


def hyperion(
    path: Path, image_models_kwargs: dict | None = None, imread_kwargs: dict | None = None
) -> SpatialData:
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs
    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 4096, 4096)
    imread_kwargs = {} if imread_kwargs is None else imread_kwargs

    files = [file for file in Path(path).iterdir() if file.suffix == ".tiff"]

    names = _get_channel_names_hyperion(files)
    image = da.concatenate(
        [imread(file, **imread_kwargs) for file in files],
        axis=0,
    )
    image = (image / image.max().compute() * 255).astype(np.uint8)

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
