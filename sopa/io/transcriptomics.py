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
from typing import Any, Optional

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray
from dask.dataframe import read_parquet
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata._logging import logger
from spatialdata.models import Image2DModel, PointsModel
from spatialdata.transformations import Affine, Identity, Scale
from spatialdata_io._constants._constants import CosmxKeys, MerscopeKeys, XeniumKeys

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
        dim=0,
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


def cosmx(
    path: str | Path,
    dataset_id: Optional[str] = None,
    fov: int | str | None = None,
    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> SpatialData:
    """
    Read *Cosmx Nanostring* data. The fields of view are stitched together, except if `fov_id` is provided.

    This function reads the following files:

        - ``<dataset_id>_`{cx.FOV_SUFFIX!r}```: Field of view file.
        - ``{cx.IMAGES_DIR!r}``: Directory containing the images.
        - ``{cx.LABELS_DIR!r}``: Directory containing the labels.

    Parameters
    ----------
    path
        Path to the root directory containing *Nanostring* files.
    dataset_id
        Name of the dataset.
    fov
        Name or number of one single field of view to be read. If a string is provided, an example of correct syntax is "F008". By default, reads all FOVs.
    imread_kwargs
        Keyword arguments passed to :func:`dask_image.imread.imread`.
    image_models_kwargs
        Keyword arguments passed to :class:`spatialdata.models.Image2DModel`.

    Returns
    -------
    :class:`spatialdata.SpatialData`
    """
    path = Path(path)

    dataset_id = _infer_dataset_id(path, dataset_id)
    fov_id, fov = _check_fov_id(fov)
    fov_locs = _read_cosmx_fov_locs(path / f"{dataset_id}_{CosmxKeys.FOV_SUFFIX}")

    ### Read image(s)
    images_dir = path / "Morphology2D"
    assert images_dir.exists(), f"Images directory not found: {images_dir}."

    if fov is None:
        image = _read_stitched_image(images_dir, fov_locs, **imread_kwargs)
        image_name = "stitched_image"
    else:
        pattern = f"*{fov_id}.TIF"
        fov_files = list(images_dir.rglob(pattern))

        assert len(fov_files), f"No file matches the pattern {pattern} inside {images_dir}"
        assert (
            len(fov_files) == 1
        ), f"Multiple files match the pattern {pattern}: {', '.join(fov_files)}"

        image = imread(images_dir / fov_files[0], **imread_kwargs)
        image_name = f"{fov}_image"

    c_coords = _cosmx_channel_names(path, image.shape[0])

    parsed_image = Image2DModel.parse(
        image, dims=("c", "y", "x"), c_coords=c_coords, **image_models_kwargs
    )

    ### Read transcripts
    transcripts_data = _read_cosmx_csv(path, dataset_id)

    if fov is None:
        transcripts_data["x"] = transcripts_data["x_global_px"] - fov_locs["xmin"].min()
        transcripts_data["y"] = transcripts_data["y_global_px"] - fov_locs["ymin"].min()
        coordinates = None
        points_name = "points"
    else:
        transcripts_data = transcripts_data[transcripts_data["fov"] == fov]
        coordinates = {"x": "x_local_px", "y": "y_local_px"}
        points_name = f"{fov}_points"

    transcripts = PointsModel.parse(
        transcripts_data,
        coordinates=coordinates,
        feature_key=CosmxKeys.TARGET_OF_TRANSCRIPT,
    )

    return SpatialData(images={image_name: parsed_image}, points={points_name: transcripts})


def _infer_dataset_id(path: Path, dataset_id: str | None) -> str:
    if isinstance(dataset_id, str):
        return dataset_id

    counts_files = list(path.rglob(f"*{CosmxKeys.TRANSCRIPTS_SUFFIX}"))
    if len(counts_files) == 1:
        found = re.match(rf"(.*)_{CosmxKeys.TRANSCRIPTS_SUFFIX}", str(counts_files[0]))
        if found:
            return found.group(1)

    raise ValueError(
        "Could not infer `dataset_id` from the name of the transcript file. Please specify it manually."
    )


def _read_cosmx_fov_locs(fov_file: Path) -> pd.DataFrame:
    assert fov_file.exists(), f"Missing field of view file: {fov_file}."

    fov_locs = pd.read_csv(fov_file, index_col=1)

    pixel_size = 0.120280945  # size of a pixel in microns

    fov_locs["xmin"] = fov_locs["X_mm"] * 1e3 / pixel_size
    fov_locs["xmax"] = 0  # will be filled when reading the images

    fov_locs["ymin"] = 0  # will be filled when reading the images
    fov_locs["ymax"] = fov_locs["Y_mm"] * 1e3 / pixel_size

    return fov_locs


def _read_stitched_image(images_dir: Path, fov_locs: pd.DataFrame, **imread_kwargs) -> da.Array:
    fov_images = {}
    pattern = re.compile(r".*_F(\d+)")
    for image_path in images_dir.iterdir():
        if image_path.suffix == ".TIF":
            fov = int(pattern.findall(image_path.name)[0])

            image = imread(image_path, **imread_kwargs)
            fov_images[fov] = da.flip(image, axis=1)

            fov_locs.loc[fov, "xmax"] = fov_locs.loc[fov, "xmin"] + image.shape[2]
            fov_locs.loc[fov, "ymin"] = fov_locs.loc[fov, "ymax"] - image.shape[1]

    for dim in ["x", "y"]:
        shift = fov_locs[f"{dim}min"].min()
        fov_locs[f"{dim}0"] = (fov_locs[f"{dim}min"] - shift).round().astype(int)
        fov_locs[f"{dim}1"] = (fov_locs[f"{dim}max"] - shift).round().astype(int)

    stitched_image = da.zeros(
        shape=(image.shape[0], fov_locs["y1"].max(), fov_locs["x1"].max()), dtype=image.dtype
    )

    for fov, im in fov_images.items():
        xmin, xmax = fov_locs.loc[fov, "x0"], fov_locs.loc[fov, "x1"]
        ymin, ymax = fov_locs.loc[fov, "y0"], fov_locs.loc[fov, "y1"]
        stitched_image[:, ymin:ymax, xmin:xmax] = im

    return stitched_image.rechunk((1, 1024, 1024))


def _check_fov_id(fov: str | int | None) -> tuple[str, int]:
    if fov is None:
        return None, None

    if isinstance(fov, int):
        return f"F{fov:0>3}", fov

    assert (
        fov[0] == "F" and len(fov) == 4 and all(c.isdigit() for c in fov[1:])
    ), f"'fov' needs to start with a F followed by three digits. Found '{fov}'."

    return fov, int(fov[1:])


def _read_cosmx_csv(path: Path, dataset_id: str) -> pd.DataFrame:
    transcripts_file = path / f"{dataset_id}_tx_file.csv.gz"

    if transcripts_file.exists():
        return pd.read_csv(transcripts_file, compression="gzip")

    transcripts_file = path / f"{dataset_id}_tx_file.csv"

    assert transcripts_file.exists(), f"Transcript file {transcripts_file} not found."

    return pd.read_csv(transcripts_file)


def _cosmx_channel_names(path: Path, n_channels: int) -> list[str]:
    channel_ids_path = path / "Morphology_ChannelID_Dictionary.txt"

    if channel_ids_path.exists():
        channel_names = list(pd.read_csv(channel_ids_path, delimiter="\t")["BiologicalTarget"])
    else:
        channel_names = [str(i) for i in range(n_channels)]
        log.warn(f"Channel file not found at {channel_ids_path}. Using {channel_names=} instead.")

    return channel_names
