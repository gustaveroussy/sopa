# Updated from spatialdata-io: https://spatialdata.scverse.org/projects/io/en/latest/
# In the future, we will completely rely on spatialdata-io

from __future__ import annotations

import logging
import re
import warnings
from pathlib import Path
from typing import Callable

import dask.array as da
import dask.dataframe as dd
import numpy as np
import xarray
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata._logging import logger
from spatialdata.models import Image2DModel, PointsModel
from spatialdata.transformations import Affine, Identity
from spatialdata_io._constants._constants import MerscopeKeys

from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


SUPPORTED_BACKENDS = ["dask_image", "rioxarray"]


def merscope(
    path: str | Path,
    backend: str = None,
    z_layers: int | list[int] | None = 3,
    region_name: str | None = None,
    slide_name: str | None = None,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
) -> SpatialData:
    """Read MERSCOPE data as a `SpatialData` object.

    This function reads the following files:
        - `detected_transcripts.csv`: transcripts locations and names
        - all the images under the `images` directory
        - `images/micron_to_mosaic_pixel_transform.csv`: affine transformation

    Args:
        path: Path to the MERSCOPE directory containing all the experiment files
        backend: Either `"dask_image"` or `"rioxarray"` (the latter uses less RAM). By default, uses `"rioxarray"` if and only if the `rioxarray` library is installed.
        z_layers: Indices of the z-layers to consider. Either one `int` index, or a list of `int` indices. If `None`, then no image is loaded. By default, only the middle layer is considered (that is, layer 3).
        region_name: Name of the region of interest, e.g., `'region_0'`. If `None` then the name of the `path` directory is used.
        slide_name: Name of the slide/run. If `None` then the name of the parent directory of `path` is used (whose name starts with a date).
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object representing the MERSCOPE experiment
    """
    assert (
        backend is None or backend in SUPPORTED_BACKENDS
    ), f"Backend '{backend} not supported. Should be one of: {', '.join(SUPPORTED_BACKENDS)}"

    path = Path(path).absolute()
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    images_dir = path / MerscopeKeys.IMAGES_DIR

    microns_to_pixels = Affine(
        np.genfromtxt(images_dir / MerscopeKeys.TRANSFORMATION_FILE),
        input_axes=("x", "y"),
        output_axes=("x", "y"),
    )

    vizgen_region = path.name if region_name is None else region_name
    slide_name = path.parent.name if slide_name is None else slide_name
    dataset_id = f"{slide_name}_{vizgen_region}"

    # Images
    images = {}

    z_layers = [z_layers] if isinstance(z_layers, int) else z_layers or []

    stainings = _get_channel_names(images_dir)
    image_transformations = {"microns": microns_to_pixels.inverse()}

    reader = _get_reader(backend)

    if stainings:
        for z_layer in z_layers:
            images[f"{dataset_id}_z{z_layer}"] = reader(
                images_dir,
                stainings,
                z_layer,
                image_models_kwargs,
                image_transformations,
                **imread_kwargs,
            )

    # Transcripts
    points = {}
    transcript_path = path / MerscopeKeys.TRANSCRIPTS_FILE
    if transcript_path.exists():
        points[f"{dataset_id}_transcripts"] = _get_points(transcript_path)
    else:
        logger.warning(
            f"Transcript file {transcript_path} does not exist. Transcripts are not loaded."
        )

    return SpatialData(points=points, images=images)


def _get_reader(backend: str | None) -> Callable:
    if backend is not None:
        return _rioxarray_load_merscope if backend == "rioxarray" else _dask_image_load_merscope
    try:
        import rioxarray  # noqa: F401

        return _rioxarray_load_merscope
    except:
        return _dask_image_load_merscope


def _get_channel_names(images_dir: Path) -> list[str]:
    exp = r"mosaic_(?P<stain>[\w|-]+[0-9]?)_z(?P<z>[0-9]+).tif"
    matches = [re.search(exp, file.name) for file in images_dir.iterdir()]

    stainings = {match.group("stain") for match in matches if match}

    return list(stainings)


def _rioxarray_load_merscope(
    images_dir: Path,
    stainings: list[str],
    z_layer: int,
    image_models_kwargs: dict,
    transformations: dict,
    **kwargs,
):
    log.info("Using rioxarray backend.")

    import rioxarray
    from rasterio.errors import NotGeoreferencedWarning

    warnings.simplefilter("ignore", category=NotGeoreferencedWarning)

    im = xarray.concat(
        [
            rioxarray.open_rasterio(
                images_dir / f"mosaic_{stain}_z{z_layer}.tif",
                chunks=image_models_kwargs["chunks"],
                **kwargs,
            )
            .rename({"band": "c"})
            .reset_coords("spatial_ref", drop=True)
            for stain in stainings
        ],
        dim="c",
    )

    return Image2DModel.parse(
        im, transformations=transformations, c_coords=stainings, **image_models_kwargs
    )


def _dask_image_load_merscope(
    images_dir: Path,
    stainings: list[str],
    z_layer: int,
    image_models_kwargs: dict,
    transformations: dict,
    **kwargs,
):
    im = da.stack(
        [
            imread(images_dir / f"mosaic_{stain}_z{z_layer}.tif", **kwargs).squeeze()
            for stain in stainings
        ],
        axis=0,
    )

    return Image2DModel.parse(
        im,
        dims=("c", "y", "x"),
        transformations=transformations,
        c_coords=stainings,
        **image_models_kwargs,
    )


def _get_points(transcript_path: Path):
    transcript_df = dd.read_csv(transcript_path)
    transcripts = PointsModel.parse(
        transcript_df,
        coordinates={"x": MerscopeKeys.GLOBAL_X, "y": MerscopeKeys.GLOBAL_Y},
        feature_key="gene",
        transformations={"microns": Identity()},
    )
    transcripts["gene"] = transcripts["gene"].astype("category")
    return transcripts
