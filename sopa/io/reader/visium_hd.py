import logging
from pathlib import Path

from spatialdata import SpatialData

from ..._constants import SopaAttrs
from ...segmentation import shapes_bounding_box
from ...utils import ensure_string_channel_names
from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def visium_hd(
    path: str | Path,
    fullres_image_file: str | Path | None = None,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
    var_names_make_unique: bool = True,
    **kwargs: int,
) -> SpatialData:
    """Read Visium HD data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.visium_hd.html).

    !!! info
        If your `fullres_image_file` is not in the `microscope_image` directory, you can specify the path to the image file using the `fullres_image_file` argument.

    Args:
        path: Path to the Visium HD directory containing all the experiment files
        fullres_image_file: Path to the full-resolution image. By default the image is searched in the `'microscope_image'` directory.
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.
        var_names_make_unique: If True, ensure that the var names are unique.
        kwargs: Additional keyword arguments passed to `spatialdata_io.visium_hd`.

    Returns:
        A `SpatialData` object representing the Xenium experiment
    """
    from spatialdata_io.readers.visium_hd import visium_hd as visium_hd_spatialdata_io

    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    del image_models_kwargs["scale_factors"]  # already set in the spatialdata_io reader

    sdata: SpatialData = visium_hd_spatialdata_io(
        path,
        fullres_image_file=fullres_image_file,
        image_models_kwargs=image_models_kwargs,
        imread_kwargs=imread_kwargs,
        **kwargs,
    )

    ensure_string_channel_names(sdata)  # Ensure that channel names are strings

    ### Add Sopa attributes to detect the spatial elements
    for key in sdata.images:
        if key.endswith("_full_image"):
            sdata.attrs[SopaAttrs.CELL_SEGMENTATION] = key
        elif key.endswith("_hires_image"):
            sdata.attrs[SopaAttrs.TISSUE_SEGMENTATION] = key

    if SopaAttrs.CELL_SEGMENTATION not in sdata.attrs:
        log.warning(
            "The full-resolution image was not found, you'll not be able to run cell segmentation.\nPlease set the `fullres_image_file` argument to the path of the full-resolution image."
        )

    for key, adata in sdata.tables.items():
        if key.endswith("_002um"):
            sdata.attrs[SopaAttrs.BINS_TABLE] = key
        if var_names_make_unique:
            adata.var_names_make_unique()

    for key in sdata.shapes:
        if key.endswith("_002um"):
            shapes_bounding_box(sdata, key)
            break

    return sdata
