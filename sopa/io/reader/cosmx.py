from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional

import dask.array as da
import pandas as pd
import tifffile
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel
from spatialdata_io._constants._constants import CosmxKeys

from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def cosmx(
    path: str | Path,
    dataset_id: Optional[str] = None,
    fov: int | str | None = None,
    read_proteins: bool = False,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
) -> SpatialData:
    """
    Read *Cosmx Nanostring* data. The fields of view are stitched together, except if `fov` is provided.

    This function reads the following files:
        - `*_fov_positions_file.csv` or `*_fov_positions_file.csv.gz`: FOV locations
        - `Morphology2D` directory: all the FOVs morphology images
        - `Morphology_ChannelID_Dictionary.txt`: Morphology channels names
        - `*_tx_file.csv.gz` or `*_tx_file.csv`: Transcripts location and names
        - If `read_proteins` is `True`, all the images under the nested `ProteinImages` directories will be read

    Args:
        path: Path to the root directory containing *Nanostring* files.
        dataset_id: Optional name of the dataset (needs to be provided if not infered).
        fov: Name or number of one single field of view to be read. If a string is provided, an example of correct syntax is "F008". By default, reads all FOVs.
        read_proteins: Whether to read the proteins or the transcripts.
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object representing the CosMX experiment
    """
    path = Path(path)
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    dataset_id = _infer_dataset_id(path, dataset_id)
    fov_locs = _read_cosmx_fov_locs(path, dataset_id)
    fov_id, fov = _check_fov_id(fov)

    protein_dir_dict = {}
    if read_proteins:
        protein_dir_dict = {
            int(protein_dir.parent.name[3:]): protein_dir
            for protein_dir in list(path.rglob("**/FOV*/ProteinImages"))
        }
        assert len(protein_dir_dict), f"No directory called 'ProteinImages' was found under {path}"

    ### Read image(s)
    images_dir = _find_dir(path, "Morphology2D")
    if fov is None:
        image, protein_names = _read_stitched_image(
            images_dir,
            fov_locs,
            protein_dir_dict,
            **imread_kwargs,
        )
        image = image.rechunk(image_models_kwargs["chunks"])
        image_name = "stitched_image"
    else:
        pattern = f"*{fov_id}.TIF"
        fov_files = list(images_dir.rglob(pattern))

        assert len(fov_files), f"No file matches the pattern {pattern} inside {images_dir}"
        assert (
            len(fov_files) == 1
        ), f"Multiple files match the pattern {pattern}: {', '.join(fov_files)}"

        image, protein_names = _read_fov_image(
            images_dir / fov_files[0], protein_dir_dict.get(fov), **imread_kwargs
        )
        image_name = f"{fov}_image"

    c_coords = _cosmx_c_coords(path, image.shape[0], protein_names)

    parsed_image = Image2DModel.parse(
        image, dims=("c", "y", "x"), c_coords=c_coords, **image_models_kwargs
    )

    if read_proteins:
        return SpatialData(images={image_name: parsed_image})

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


def _read_fov_image(
    morphology_path: Path, protein_path: Path | None, **imread_kwargs
) -> tuple[da.Array, list[str] | None]:
    image = imread(morphology_path, **imread_kwargs)

    protein_names = None
    if protein_path is not None:
        protein_image, protein_names = _read_protein_fov(protein_path)
        image = da.concatenate([image, protein_image], axis=0)

    return image, protein_names


def _infer_dataset_id(path: Path, dataset_id: str | None) -> str:
    if isinstance(dataset_id, str):
        return dataset_id

    for suffix in [".csv", ".csv.gz"]:
        counts_files = list(path.rglob(f"*_fov_positions_file{suffix}"))

        if len(counts_files) == 1:
            found = re.match(rf"(.*)_fov_positions_file{suffix}", str(counts_files[0]))
            if found:
                return found.group(1)

    raise ValueError(
        "Could not infer `dataset_id` from the name of the transcript file. Please specify it manually."
    )


def _read_cosmx_fov_locs(path: Path, dataset_id: str) -> pd.DataFrame:
    fov_file = path / f"{dataset_id}_fov_positions_file.csv"

    if not fov_file.exists():
        fov_file = path / f"{dataset_id}_fov_positions_file.csv.gz"

    assert fov_file.exists(), f"Missing field of view file: {fov_file}"

    fov_locs = pd.read_csv(fov_file, index_col=1)

    pixel_size = 0.120280945  # size of a pixel in microns

    fov_locs["xmin"] = fov_locs["X_mm"] * 1e3 / pixel_size
    fov_locs["xmax"] = 0  # will be filled when reading the images

    fov_locs["ymin"] = 0  # will be filled when reading the images
    fov_locs["ymax"] = fov_locs["Y_mm"] * 1e3 / pixel_size

    return fov_locs


def _read_stitched_image(
    images_dir: Path, fov_locs: pd.DataFrame, protein_dir_dict: dict, **imread_kwargs
) -> tuple[da.Array, list[str] | None]:
    log.warn("Image stitching is currently experimental")

    fov_images = {}
    protein_names = None
    pattern = re.compile(r".*_F(\d+)")
    for image_path in images_dir.iterdir():
        if image_path.suffix == ".TIF":
            fov = int(pattern.findall(image_path.name)[0])

            image, protein_names = _read_fov_image(
                image_path, protein_dir_dict.get(fov), **imread_kwargs
            )

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

    return stitched_image, protein_names


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


def _cosmx_c_coords(path: Path, n_channels: int, protein_names: list[str] | None) -> list[str]:
    channel_ids_path = path / "Morphology_ChannelID_Dictionary.txt"

    if channel_ids_path.exists():
        channel_names = list(pd.read_csv(channel_ids_path, delimiter="\t")["BiologicalTarget"])
    else:
        n_channels = n_channels - len(protein_names) if protein_names is not None else n_channels
        channel_names = [str(i) for i in range(n_channels)]
        log.warn(f"Channel file not found at {channel_ids_path}, using {channel_names=} instead.")

    if protein_names is not None:
        channel_names += protein_names

    return channel_names


def _find_dir(path: Path, name: str):
    if (path / name).is_dir():
        return path / name

    paths = list(path.rglob(f"**/{name}"))
    assert len(paths) == 1, f"Found {len(paths)} path(s) with name {name} inside {path}"

    return paths[0]


def _get_cosmx_protein_name(image_path: Path) -> str:
    with tifffile.TiffFile(image_path) as tif:
        description = tif.pages[0].description
        substrings = re.findall(r'"DisplayName": "(.*?)",', description)
        return substrings[0]


def _read_protein_fov(protein_dir: Path) -> tuple[da.Array, list[str]]:
    images_paths = list(protein_dir.rglob("*.TIF"))
    protein_image = da.concatenate([imread(image_path) for image_path in images_paths], axis=0)
    channel_names = [_get_cosmx_protein_name(image_path) for image_path in images_paths]

    return protein_image, channel_names
