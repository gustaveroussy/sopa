import logging
from pathlib import Path
from typing import Literal

from spatialdata import SpatialData

from ..._constants import SopaAttrs
from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def merscope(
    path: str | Path,
    backend: Literal["dask_image", "rioxarray"] | None = None,
    z_layers: int | list[int] | None = 3,
    region_name: str | None = None,
    slide_name: str | None = None,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
    **kwargs: int,
) -> SpatialData:
    """Read MERSCOPE data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.merscope.html).

    This function reads the following files:
        - `detected_transcripts.csv`: transcripts locations and names
        - all the images under the `images` directory
        - `images/micron_to_mosaic_pixel_transform.csv`: affine transformation

    Args:
        path: Path to the MERSCOPE directory containing all the experiment files
        backend: Either `"dask_image"` or `"rioxarray"` (the latter uses less RAM, but requires `rioxarray` to be installed). By default, uses `"rioxarray"` if and only if the `rioxarray` library is installed.
        z_layers: Indices of the z-layers to consider. Either one `int` index, or a list of `int` indices. If `None`, then no image is loaded. By default, only the middle layer is considered (that is, layer 3).
        region_name: Name of the region of interest, e.g., `'region_0'`. If `None` then the name of the `path` directory is used.
        slide_name: Name of the slide/run. If `None` then the name of the parent directory of `path` is used (whose name starts with a date).
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object representing the MERSCOPE experiment
    """
    from spatialdata_io.readers.merscope import merscope as merscope_spatialdata_io

    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    sdata = merscope_spatialdata_io(
        path,
        backend=backend,
        z_layers=z_layers,
        region_name=region_name,
        slide_name=slide_name,
        image_models_kwargs=image_models_kwargs,
        imread_kwargs=imread_kwargs,
        cells_boundaries=False,
        cells_table=False,
        **kwargs,
    )

    ### Add Sopa attributes to detect the spatial elements
    if z_layers is not None:
        if not isinstance(z_layers, int) and len(z_layers) == 1:
            z_layers = z_layers[0]
        if isinstance(z_layers, int):
            for key in sdata.images.keys():
                if key.endswith(f"_z{z_layers}"):
                    sdata.attrs[SopaAttrs.CELL_SEGMENTATION] = key
        else:
            log.warning(
                f"Multiple z-layers provided: {z_layers}. Not deciding which image should be used for cell segmentation."
            )

    for key in sdata.points.keys():
        if key.endswith("_transcripts"):
            sdata.attrs[SopaAttrs.TRANSCRIPTS] = key

    return sdata
