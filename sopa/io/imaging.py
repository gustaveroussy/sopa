import logging
import re
from pathlib import Path

import dask.array as da
import tifffile as tf
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

log = logging.getLogger(__name__)


def get_channel_name(description):
    return re.search(r"<Name>(.*?)</Name>", description).group(1)


def qptiff(
    path: Path, channels_renaming: dict | None = None, image_models_kwargs: dict | None = None
) -> SpatialData:
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs
    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 4096, 4096)

    with tf.TiffFile(path) as tif:
        page_series = tif.series[0]
        names = [get_channel_name(page.description) for page in page_series]

        log.info(f"Found channel names {names}")

        if len(channels_renaming):
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
