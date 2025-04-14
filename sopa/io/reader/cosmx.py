import logging
import re
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import tifffile
import xarray as xr
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel

from ..._constants import SopaAttrs
from .utils import _deduplicate_names, _default_image_kwargs

log = logging.getLogger(__name__)


def cosmx(
    path: str | Path,
    dataset_id: str | None = None,
    fov: int | None = None,
    read_proteins: bool = False,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
    flip_image: bool | None = None,
) -> SpatialData:
    """
    Read *Cosmx Nanostring* data. The fields of view are stitched together, except if `fov` is provided.

    This function reads the following files:
        - `*_fov_positions_file.csv` or `*_fov_positions_file.csv.gz`: FOV locations
        - `Morphology2D` directory: all the FOVs morphology images
        - `*_tx_file.csv.gz` or `*_tx_file.csv`: Transcripts location and names
        - If `read_proteins` is `True`, all the images under the nested `ProteinImages` directories will be read

        These files must be exported as flat files in AtomX. That is: within a study, click on "Export" and then select files from the "Flat CSV Files" section (transcripts flat and FOV position flat).

    Args:
        path: Path to the root directory containing *Nanostring* files.
        dataset_id: Optional name of the dataset (needs to be provided if not inferred).
        fov: Number of one single field of view to be read. If not provided, reads all FOVs and create a stitched image.
        read_proteins: Whether to read the proteins or the transcripts.
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.
        flip_image: For some buggy AtomX exports, `flip_image=True` has to be used for stitching. By default, the value is inferred based on the transcript file. See [this](https://github.com/gustaveroussy/sopa/issues/231) issue.

    Returns:
        A `SpatialData` object representing the CosMX experiment
    """
    path = Path(path)
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    dataset_id = _infer_dataset_id(path, dataset_id)
    fov_locs = _read_fov_locs(path, dataset_id)

    protein_dir_dict = {}
    if read_proteins:
        protein_dir_dict = {
            int(protein_dir.parent.name[3:]): protein_dir for protein_dir in list(path.rglob("**/FOV*/ProteinImages"))
        }
        assert len(protein_dir_dict), f"No directory called 'ProteinImages' was found under {path}"

    if flip_image is None and not read_proteins:
        flip_image = _infer_flip_image(path, dataset_id)

    ### Read image(s)
    images_dir = _find_dir(path, "Morphology2D")
    morphology_coords = _cosmx_morphology_coords(images_dir)

    if fov is None:
        image, c_coords = _read_stitched_image(
            images_dir,
            fov_locs,
            protein_dir_dict,
            morphology_coords,
            flip_image,
            **imread_kwargs,
        )
        image_name = "stitched_image"
    else:
        log.info(f"Reading single FOV ({fov}), the image will not be stitched")
        fov_file = _find_matching_fov_file(images_dir, fov)

        image, c_coords = _read_fov_image(fov_file, protein_dir_dict.get(fov), morphology_coords, **imread_kwargs)
        image_name = f"F{fov:0>5}_image"

    parsed_image = Image2DModel.parse(image, dims=("c", "y", "x"), c_coords=c_coords, **image_models_kwargs)

    if read_proteins:
        return SpatialData(images={image_name: parsed_image}, attrs={SopaAttrs.CELL_SEGMENTATION: image_name})

    ### Read transcripts
    transcripts_data = _read_transcripts_csv(path, dataset_id)

    if fov is None:
        transcripts_data["x"] = transcripts_data["x_global_px"] - fov_locs["xmin"].min()
        transcripts_data["y"] = transcripts_data["y_global_px"] - fov_locs["ymin"].min()
        coordinates = None
        points_name = "points"
    else:
        transcripts_data = transcripts_data[transcripts_data["fov"] == fov]
        coordinates = {"x": "x_local_px", "y": "y_local_px"}
        points_name = f"F{fov:0>5}_points"

    from spatialdata_io._constants._constants import CosmxKeys

    transcripts = PointsModel.parse(
        transcripts_data,
        coordinates=coordinates,
        feature_key=CosmxKeys.TARGET_OF_TRANSCRIPT,
    )

    return SpatialData(
        images={image_name: parsed_image},
        points={points_name: transcripts},
        attrs={
            SopaAttrs.CELL_SEGMENTATION: image_name,
            SopaAttrs.TRANSCRIPTS: points_name,
            SopaAttrs.PRIOR_TUPLE_KEY: ("unique_cell_id", 0),
        },
    )


def _infer_dataset_id(path: Path, dataset_id: str | None) -> str:
    if isinstance(dataset_id, str):
        return dataset_id

    for suffix in [".csv", ".csv.gz"]:
        counts_files = list(path.rglob(f"[!\\.]*_fov_positions_file{suffix}"))

        if len(counts_files) == 1:
            found = re.match(rf"(.*)_fov_positions_file{suffix}", counts_files[0].name)
            if found:
                return found.group(1)

    raise ValueError("Could not infer `dataset_id` from the name of the transcript file. Please specify it manually.")


def _read_fov_image(
    morphology_path: Path, protein_path: Path | None, morphology_coords: list[str], **imread_kwargs
) -> tuple[da.Array, list[str]]:
    image = imread(morphology_path, **imread_kwargs)

    protein_names = []
    if protein_path is not None:
        protein_image, protein_names = _read_protein_fov(protein_path)
        image = da.concatenate([image, protein_image], axis=0)

    return image, _deduplicate_names(morphology_coords + protein_names)


def _read_fov_locs(path: Path, dataset_id: str) -> pd.DataFrame:
    fov_file = path / f"{dataset_id}_fov_positions_file.csv"

    if not fov_file.exists():
        fov_file = path / f"{dataset_id}_fov_positions_file.csv.gz"

    assert fov_file.exists(), f"Missing field of view file: {fov_file}"

    fov_locs = pd.read_csv(fov_file)
    fov_locs["xmax"] = 0.0  # will be filled when reading the images
    fov_locs["ymax"] = 0.0  # will be filled when reading the images

    fov_key, x_key, y_key, scale_factor = "fov", "x_global_px", "y_global_px", 1

    if not np.isin([fov_key, x_key, y_key], fov_locs.columns).all():  # try different column names
        fov_key, x_key, y_key = "FOV", "X_mm", "Y_mm"
        scale_factor = 1e3 / 0.120280945  # CosMX milimeters to pixels

        assert np.isin([fov_key, x_key, y_key], fov_locs.columns).all(), (
            f"The file {fov_file} must contain the following columns: {fov_key}, {x_key}, {y_key}. Consider using a different export module."
        )

    fov_locs.index = fov_locs[fov_key]
    fov_locs["xmin"] = fov_locs[x_key] * scale_factor
    fov_locs["ymin"] = fov_locs[y_key] * scale_factor

    return fov_locs


def _read_stitched_image(
    images_dir: Path,
    fov_locs: pd.DataFrame,
    protein_dir_dict: dict,
    morphology_coords: list[str],
    flip_image: int,
    **imread_kwargs,
) -> tuple[da.Array, list[str] | None]:
    fov_images = {}
    c_coords_dict = {}
    pattern = re.compile(r".*_F(\d+)")
    for image_path in images_dir.iterdir():
        if image_path.suffix == ".TIF":
            fov = int(pattern.findall(image_path.name)[0])

            image, c_coords = _read_fov_image(image_path, protein_dir_dict.get(fov), morphology_coords, **imread_kwargs)

            c_coords_dict[fov] = c_coords

            fov_images[fov] = da.flip(image, axis=1)

            fov_locs.loc[fov, "xmax"] = fov_locs.loc[fov, "xmin"] + image.shape[2]

            if flip_image:
                fov_locs.loc[fov, "ymax"] = fov_locs.loc[fov, "ymin"]
                fov_locs.loc[fov, "ymin"] = fov_locs.loc[fov, "ymax"] - image.shape[1]
            else:
                fov_locs.loc[fov, "ymax"] = fov_locs.loc[fov, "ymin"] + image.shape[1]

    for dim in ["x", "y"]:
        shift = fov_locs[f"{dim}min"].min()
        fov_locs[f"{dim}0"] = (fov_locs[f"{dim}min"] - shift).round().astype(int)
        fov_locs[f"{dim}1"] = (fov_locs[f"{dim}max"] - shift).round().astype(int)

    c_coords = list(set.union(*[set(names) for names in c_coords_dict.values()]))

    height, width = fov_locs["y1"].max(), fov_locs["x1"].max()
    stitched_image = da.zeros(shape=(len(c_coords), height, width), dtype=image.dtype)
    stitched_image = xr.DataArray(stitched_image, dims=("c", "y", "x"), coords={"c": c_coords})

    for fov, im in fov_images.items():
        xmin, xmax = fov_locs.loc[fov, "x0"], fov_locs.loc[fov, "x1"]
        ymin, ymax = fov_locs.loc[fov, "y0"], fov_locs.loc[fov, "y1"]

        if flip_image:
            y_slice, x_slice = slice(height - ymax, height - ymin), slice(width - xmax, width - xmin)
        else:
            y_slice, x_slice = slice(ymin, ymax), slice(xmin, xmax)

        stitched_image.loc[{"c": c_coords_dict[fov], "y": y_slice, "x": x_slice}] = im

        if len(c_coords_dict[fov]) < len(c_coords):
            log.warning(f"Missing channels ({len(c_coords) - len(c_coords_dict[fov])}) for FOV {fov}")

    return stitched_image.data, c_coords


def _find_matching_fov_file(images_dir: Path, fov: str | int) -> Path:
    assert isinstance(fov, int), "Expected `fov` to be an integer"

    pattern = re.compile(rf".*_F0*{fov}\.TIF")
    fov_files = [file for file in images_dir.rglob("*") if pattern.match(file.name)]

    assert len(fov_files), f"No file matches the pattern {pattern} inside {images_dir}"
    assert len(fov_files) == 1, f"Multiple files match the pattern {pattern}: {', '.join(map(str, fov_files))}"

    return fov_files[0]


def _read_transcripts_csv(path: Path, dataset_id: str, nrows: int | None = None) -> pd.DataFrame:
    transcripts_file = path / f"{dataset_id}_tx_file.csv.gz"

    if transcripts_file.exists():
        df = pd.read_csv(transcripts_file, compression="gzip", nrows=nrows)
    else:
        transcripts_file = path / f"{dataset_id}_tx_file.csv"
        assert transcripts_file.exists(), f"Transcript file {transcripts_file} not found."
        df = pd.read_csv(transcripts_file, nrows=nrows)

    TRANSCRIPT_COLUMNS = ["x_global_px", "y_global_px", "target"]
    assert np.isin(TRANSCRIPT_COLUMNS, df.columns).all(), (
        f"The file {transcripts_file} must contain the following columns: {', '.join(TRANSCRIPT_COLUMNS)}. Consider using a different export module."
    )

    df["unique_cell_id"] = df["fov"] * (df["cell_ID"].max() + 1) * (df["cell_ID"] > 0) + df["cell_ID"]

    return df


def _find_dir(path: Path, name: str):
    if (path / name).is_dir():
        return path / name

    paths = list(path.rglob(f"**/{name}"))
    assert len(paths) == 1, f"Found {len(paths)} path(s) with name {name} inside {path}"

    return paths[0]


def _cosmx_morphology_coords(images_dir: Path) -> list[str]:
    images_paths = list(images_dir.glob("*.TIF"))
    assert len(images_paths) > 0, f"Expected to find images inside {images_dir}"

    with tifffile.TiffFile(images_paths[0]) as tif:
        description = tif.pages[0].description

        substrings = re.findall(r'"BiologicalTarget": "(.*?)",', description)
        channels = re.findall(r'"ChannelId": "(.*?)",', description)
        channel_order = list(re.findall(r'"ChannelOrder": "(.*?)",', description)[0])

        return [substrings[channels.index(x)] for x in channel_order if x in channels]


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


def _infer_flip_image(path: Path, dataset_id: str) -> bool:
    df_ = _read_transcripts_csv(path, dataset_id, nrows=100)

    fov = df_["fov"].value_counts().index[0]
    df_ = df_[df_["fov"] == fov].sort_values("y_global_px")

    assert len(df_) > 1, f"Not transcripts in {fov=} to infer `flip_image`. Please provide `flip_image` manually."

    # in recent AtomX exports, y_local_px is negatively correlated with y_global_px
    flip_image = df_["y_local_px"].iloc[0] > df_["y_local_px"].iloc[-1]

    log.info(
        f"Inferring argument {flip_image=}. If the image stitching is wrong, please add a comment to https://github.com/gustaveroussy/sopa/issues/231"
    )

    return flip_image
