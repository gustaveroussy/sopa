# Readers for spatial-transcriptomics technologies
# Updated from spatialdata-io: https://spatialdata.scverse.org/projects/io/en/latest/
# In the future, we will completely rely on spatialdata-io (when stable enough)

from __future__ import annotations

import json
import logging
import re
import warnings
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import Any

import dask.dataframe as dd
import numpy as np
import spatialdata_io
import xarray
from dask import array as da
from dask.dataframe import read_parquet
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata._logging import logger
from spatialdata.models import Image2DModel, PointsModel
from spatialdata.transformations import Affine, Identity, Scale
from spatialdata_io._constants._constants import MerscopeKeys, XeniumKeys

log = logging.getLogger(__name__)


def _get_channel_names(images_dir: Path) -> list[str]:
    exp = r"mosaic_(?P<stain>[\w|-]+[0-9]?)_z(?P<z>[0-9]+).tif"
    matches = [re.search(exp, file.name) for file in images_dir.iterdir()]

    stainings = {match.group("stain") for match in matches if match}

    return list(stainings)


SUPPORTED_BACKENDS = ["dask_image", "rioxarray"]


def merscope(
    path: str | Path,
    z_layers: int | list[int] | None = 3,
    region_name: str | None = None,
    slide_name: str | None = None,
    backend: str = "dask_image",
    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> SpatialData:
    """Read MERSCOPE data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.merscope.html).

    Args:
        path: Path to the MERSCOPE directory containing all the experiment files
        backend: Either `dask_image` or `rioxarray` (the latter should use less RAM, but it is still experimental)
        **kwargs: See link above.

    Returns:
        A `SpatialData` object representing the MERSCOPE experiment
    """
    assert (
        backend in SUPPORTED_BACKENDS
    ), f"Backend '{backend} not supported. Should be one of: {', '.join(SUPPORTED_BACKENDS)}"

    if backend == "rioxarray":
        log.info("Using experimental rioxarray backend.")

    if "chunks" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["chunks"] = (1, 1024, 1024)
    if "scale_factors" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    path = Path(path).absolute()
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
    reader = _rioxarray_load_merscope if backend == "rioxarray" else _dask_image_load_merscope

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


def _rioxarray_load_merscope(
    images_dir: Path,
    stainings: list[str],
    z_layer: int,
    image_models_kwargs: dict,
    transformations: dict,
    **kwargs,
):
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
        transformations={"microns": Identity()},
    )
    transcripts["gene"] = transcripts["gene"].astype("category")
    return transcripts


def xenium(
    path: str | Path,
    imread_kwargs=MappingProxyType({}),
    image_models_kwargs=MappingProxyType({}),
) -> SpatialData:
    """Read Xenium data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.xenium.html).

    Args:
        path: Path to the Xenium directory containing all the experiment files
        imread_kwargs: See link above.
        image_models_kwargs:See link above.

    Returns:
        A `SpatialData` object representing the Xenium experiment
    """
    if "chunks" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["chunks"] = (1, 1024, 1024)
    if "scale_factors" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    path = Path(path)
    with open(path / XeniumKeys.XENIUM_SPECS) as f:
        specs = json.load(f)

    points = {"transcripts": _get_points_xenium(path, specs)}

    images = {
        "morphology_mip": _get_images_xenium(
            path,
            XeniumKeys.MORPHOLOGY_MIP_FILE,
            imread_kwargs,
            image_models_kwargs,
        )
    }

    return SpatialData(images=images, points=points)


def _get_points_xenium(path: Path, specs: dict[str, Any]):
    table = read_parquet(path / XeniumKeys.TRANSCRIPTS_FILE)
    table["feature_name"] = table["feature_name"].apply(
        lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x),
        meta=("feature_name", "object"),
    )

    transform = Scale([1.0 / specs["pixel_size"], 1.0 / specs["pixel_size"]], axes=("x", "y"))
    points = PointsModel.parse(
        table,
        coordinates={
            "x": XeniumKeys.TRANSCRIPTS_X,
            "y": XeniumKeys.TRANSCRIPTS_Y,
            "z": XeniumKeys.TRANSCRIPTS_Z,
        },
        feature_key=XeniumKeys.FEATURE_NAME,
        instance_key=XeniumKeys.CELL_ID,
        transformations={"global": transform},
    )
    return points


def _get_images_xenium(
    path: Path,
    file: str,
    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    image = imread(path / file, **imread_kwargs)
    return Image2DModel.parse(
        image,
        transformations={"global": Identity()},
        dims=("c", "y", "x"),
        c_coords=list(map(str, range(len(image)))),
        **image_models_kwargs,
    )


def cosmx(path: str, **kwargs: int) -> SpatialData:
    """Alias to the [spatialdata-io reader](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.cosmx.html).

    Args:
        path: Path to the CosMX data directory
        **kwargs: See link above.

    Returns:
        A `SpatialData` object representing the CosMX experiment
    """
    # TODO: add stitching + set chunksize to 1024
    return spatialdata_io.cosmx(path, **kwargs)
