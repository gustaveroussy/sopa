from __future__ import annotations

import logging
from pathlib import Path

from spatialdata import SpatialData
from spatialdata_io.readers.visium_hd import visium_hd as visium_hd_spatialdata_io

from ..._constants import SopaAttrs
from ...utils import ensure_string_channel_names
from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def visium_hd(
    path: str | Path,
    bin_size: int | list[int] | None = 2,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
    **kwargs: int,
) -> SpatialData:
    """Read Visium HD data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.visium_hd.html).

    Args:
        path: Path to the Visium HD directory containing all the experiment files
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object representing the Xenium experiment
    """
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    del image_models_kwargs["scale_factors"]  # already set in the spatialdata_io reader

    sdata: SpatialData = visium_hd_spatialdata_io(
        path,
        bin_size=bin_size,
        image_models_kwargs=image_models_kwargs,
        imread_kwargs=imread_kwargs,
        **kwargs,
    )

    ensure_string_channel_names(sdata)  # Ensure that channel names are strings

    ### Add Sopa attributes to detect the spatial elements
    for key in sdata.images.keys():
        if key.endswith("_full_image"):
            sdata.attrs[SopaAttrs.CELL_SEGMENTATION] = key
        elif key.endswith("_hires_image"):
            sdata.attrs[SopaAttrs.TISSUE_SEGMENTATION] = key

    for key in sdata.tables.keys():
        if key.endswith("_002um"):
            sdata.attrs[SopaAttrs.BINS_TABLE] = key

    return sdata
