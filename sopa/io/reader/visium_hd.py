from __future__ import annotations

import logging
from pathlib import Path

from spatialdata import SpatialData
from spatialdata_io.readers.visium_hd import visium_hd as visium_hd_spatialdata_io

from ..._constants import SopaAttrs
from ..._sdata import _update_spatialdata_attrs
from ...utils import string_channel_names
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

    string_channel_names(sdata)  # Ensure that channel names are strings

    ### Add Sopa attributes to detect the spatial elements
    for key, image in sdata.images.items():
        if key.endswith("_full_image"):
            _update_spatialdata_attrs(image, {SopaAttrs.CELL_SEGMENTATION: True})
        elif key.endswith("_hires_image"):
            _update_spatialdata_attrs(image, {SopaAttrs.TISSUE_SEGMENTATION: True})

    for key, geo_df in sdata.shapes.items():
        if key.endswith("_002um"):
            _update_spatialdata_attrs(geo_df, {SopaAttrs.BINS_AGGREGATION: True})

    for key, table in sdata.tables.items():
        if key.endswith("_002um"):
            _update_spatialdata_attrs(table, {SopaAttrs.BINS_TABLE: True})

    return sdata
