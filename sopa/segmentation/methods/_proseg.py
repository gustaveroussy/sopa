import gzip
import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd
import spatialdata
from anndata import AnnData
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from ..._constants import SopaAttrs, SopaKeys
from ...aggregation.aggregation import add_standardized_table
from ...utils import delete_transcripts_patches_dirs, get_feature_key, get_transcripts_patches_dirs
from .._transcripts import _check_transcript_patches
from ._utils import _get_executable_path

log = logging.getLogger(__name__)


def proseg(
    sdata: SpatialData,
    delete_cache: bool = True,
    command_line_suffix: str = "",
    infer_presets: bool = True,
    key_added: str = SopaKeys.PROSEG_BOUNDARIES,
):
    """Run [`proseg`](https://github.com/dcjones/proseg) segmentation on a SpatialData object, and add the corresponding cell boundaries and `AnnData` table with counts.

    !!! warning "Proseg installation"
        Make sure to install [`proseg`](https://github.com/dcjones/proseg) separately before running this function.

    !!! info "Proseg usage specificities"
        Contrary to most other segmentation tools, `proseg` will only run on one patch. I.e., you need
        to run [`sopa.make_transcript_patches`](../patches/#sopa.make_transcript_patches) with `patch_width=None` and a `prior_shapes_key` before running `proseg`.

        Also, note that aggregation is not necessary after running `proseg`.

    Args:
        sdata: A `SpatialData` object.
        delete_cache: Whether to delete the cache after segmentation.
        command_line_suffix: Optional suffix to add to the proseg command line.
        infer_presets: Whether to infer the proseg presets based on the columns of the transcripts dataframe.
        key_added: Name of the shapes element to be added to `sdata.shapes`.
    """
    _check_transcript_patches(sdata)

    points_key = sdata[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.POINTS_KEY].iloc[0]

    patches_dirs = get_transcripts_patches_dirs(sdata)
    assert len(patches_dirs) == 1, (
        "Proseg is fast enough to work on a single patch. Re-run `sopa.make_transcript_patches` with `patch_width=None` and a `prior_shapes_key`."
    )
    patch_dir = Path(patches_dirs[0])

    proseg_command = _get_proseg_command(sdata, points_key, command_line_suffix, infer_presets)

    _run_proseg(proseg_command, patch_dir)
    adata, geo_df = _read_proseg(sdata, patch_dir, points_key)

    add_standardized_table(sdata, adata, geo_df, key_added, SopaKeys.TABLE)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added

    if delete_cache:
        delete_transcripts_patches_dirs(sdata)

    log.info("Proseg table and boundaries added (running `sopa.aggregate` is not mandatory).")


def _run_proseg(proseg_command: str, patch_dir: str | Path):
    import subprocess

    log.info(f"Running proseg with command: `{proseg_command}`")

    result = subprocess.run(
        proseg_command,
        cwd=patch_dir,
        shell=True,
        capture_output=False,
    )

    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            returncode=result.returncode,
            cmd=proseg_command,
            output=result.stdout,
            stderr=result.stderr,
        )


def _get_proseg_command(
    sdata: SpatialData,
    points_key: str,
    command_line_suffix: str,
    infer_presets: bool,
) -> str:
    proseg_executable = _get_executable_path("proseg", ".cargo")

    assert SopaKeys.PRIOR_SHAPES_KEY in sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES], (
        "Proseg requires a prior. Re-run `sopa.make_transcript_patches` with a `prior_shapes_key`."
    )

    prior_shapes_key = sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.PRIOR_SHAPES_KEY].iloc[0]

    feature_key = get_feature_key(sdata[points_key], raise_error=True)

    use_zarr = _use_zarr_output(proseg_executable)

    if infer_presets:
        command_line_suffix = _add_presets(command_line_suffix, sdata[points_key].columns)

    return f"{proseg_executable} transcripts.csv -x x -y y -z z --gene-column {feature_key} --cell-id-column {prior_shapes_key} --cell-id-unassigned 0 {'--exclude-spatialdata-transcripts' if use_zarr else ''} {command_line_suffix}"


def _add_presets(command_line_suffix: str, columns: list[str]) -> str:
    if "fov-column" not in command_line_suffix:
        for column in columns:
            if column in ["fov", "fov_name"]:
                command_line_suffix += f" --fov-column {column}"
                break

    if "qv" in columns and "--qv-column" not in command_line_suffix:
        command_line_suffix += " --qv-column qv"

    if "qv" in columns and "--min-qv" not in command_line_suffix:
        command_line_suffix += " --min-qv 20.0"

    if "compartment-column" not in command_line_suffix:
        if "overlaps_nucleus" in columns:
            command_line_suffix += " --compartment-column overlaps_nucleus --compartment-nuclear 1"
        elif "CellComp" in columns:
            command_line_suffix += " --compartment-column CellComp --compartment-nuclear Nuclear"

    return command_line_suffix


def _read_proseg(sdata: SpatialData, patch_dir: Path, points_key: str) -> tuple[AnnData, gpd.GeoDataFrame]:
    zarr_output = patch_dir / "proseg-output.zarr"

    if zarr_output.exists():  # in versions >= 3.0.0
        _sdata = spatialdata.read_zarr(zarr_output)
        adata, geo_df = _sdata.tables["table"], _sdata.shapes["cell_boundaries"]

        del geo_df.attrs["transform"]
    else:
        counts = pd.read_csv(patch_dir / "expected-counts.csv.gz")

        obs = pd.read_csv(patch_dir / "cell-metadata.csv.gz")
        obs.index = obs.index.map(str)

        adata = AnnData(counts, obs=obs)

        with gzip.open(patch_dir / "cell-polygons.geojson.gz", "rb") as f:
            geo_df = gpd.read_file(f)

    geo_df.crs = None

    transformations = get_transformation(sdata[points_key], get_all=True).copy()
    geo_df = ShapesModel.parse(geo_df, transformations=transformations)

    return adata, geo_df


def _use_zarr_output(proseg_executable_path: str) -> bool:
    import subprocess

    from packaging.version import InvalidVersion, Version

    result = subprocess.run(f"{proseg_executable_path} --version", shell=True, capture_output=True, text=True)

    try:
        return Version(result.stdout.split()[1]) >= Version("3.0.0")
    except InvalidVersion:
        log.warning(
            "Could not parse the version of proseg. Assumes proseg<3.0.0. We may assume fallback to a different version in the future."
        )
        return False
