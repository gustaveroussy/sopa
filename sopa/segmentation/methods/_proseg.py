import gzip
import logging
import shutil
from pathlib import Path

import geopandas as gpd
import pandas as pd
from anndata import AnnData
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from ..._constants import SopaAttrs, SopaKeys
from ...aggregation.aggregation import add_standardized_table
from ...utils import (
    delete_transcripts_patches_dirs,
    get_feature_key,
    get_transcripts_patches_dirs,
)
from .._transcripts import _check_transcript_patches

log = logging.getLogger(__name__)


def proseg(
    sdata: SpatialData,
    delete_cache: bool = True,
    command_line_suffix: str = "",
    key_added: str = SopaKeys.PROSEG_BOUNDARIES,
):
    """Run [`proseg`](https://github.com/dcjones/proseg) segmentation on a SpatialData object, and add the corresponding cell boundaries and `AnnData` table with counts.

    !!! warning "Proseg installation"
        Make sure to install [`proseg`](https://github.com/dcjones/proseg) separately before running this function.

    !!! info "Proseg usage specificities"
        Contrary to most other segmentation tools, `proseg` will only
        run on one patch. I.e., you need to run
        [`sopa.make_transcript_patches`](../patches/#sopa.make_transcript_patches)
        with `patch_width=None` and (optionally) a `prior_shapes_key` before running `proseg`.

        Also, note that aggregation is not necessary after running `proseg`.

    Args:
        sdata: A `SpatialData` object.
        delete_cache: Whether to delete the cache after segmentation.
        command_line_suffix: Optional suffix to add to the proseg command line.
        key_added: Name of the shapes element to be added to `sdata.shapes`.
    """
    _check_transcript_patches(sdata)

    points_key = sdata[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.POINTS_KEY].iloc[0]

    patches_dirs = get_transcripts_patches_dirs(sdata)
    assert (
        len(patches_dirs) == 1
    ), "Proseg is fast enough to work on a single patch. Re-run `sopa.make_transcript_patches` with `patch_width=None` and a `prior_shapes_key`."
    patch_dir = Path(patches_dirs[0])

    proseg_command = _get_proseg_command(sdata, points_key, command_line_suffix)

    _run_proseg(proseg_command, patch_dir)
    adata, geo_df = _read_proseg(sdata, patch_dir, points_key)

    add_standardized_table(sdata, adata, geo_df, key_added, SopaKeys.TABLE)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added

    if delete_cache:
        delete_transcripts_patches_dirs(sdata)

    log.info("Proseg table and boundaries added (running `sopa.aggregate` is not mandatory).")


def _run_proseg(proseg_command: str, patch_dir: str | Path):
    import subprocess

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


def _get_proseg_executable_path() -> Path | str:
    if shutil.which("proseg") is not None:
        return "proseg"

    default_path = Path.home() / ".cargo" / "bin" / "proseg"
    if default_path.exists():
        return default_path

    bin_path = Path.home() / ".local" / "bin" / "proseg"
    raise FileNotFoundError(
        f"Please install proseg and ensure that either `{default_path}` executes proseg, or that `proseg` is an existing command (add it to your PATH, or create a symlink at {bin_path})."
    )


def _get_proseg_command(sdata: SpatialData, points_key: str, command_line_suffix: str) -> str:
    proseg_executable = _get_proseg_executable_path()

    if (
        SopaKeys.PRIOR_SHAPES_KEY
        in sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES]
    ):
        prior_shapes_key = sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES][
            SopaKeys.PRIOR_SHAPES_KEY
        ].iloc[0]
        prior_command = f"--cell-id-column {prior_shapes_key}"
    else:
        prior_shapes_key = None
        prior_command = ""

    feature_key = get_feature_key(sdata[points_key], raise_error=True)

    return f"{proseg_executable} transcripts.csv -x x -y y -z z --gene-column {feature_key} --cell-id-column {prior_shapes_key} --cell-id-unassigned 0 {command_line_suffix}"


def _read_proseg(sdata: SpatialData, patch_dir: Path, points_key: str) -> tuple[AnnData, gpd.GeoDataFrame]:
    counts = pd.read_csv(patch_dir / "expected-counts.csv.gz")

    obs = pd.read_csv(patch_dir / "cell-metadata.csv.gz")
    obs.index = obs.index.map(str)

    adata = AnnData(counts, obs=obs)

    with gzip.open(patch_dir / "cell-polygons.geojson.gz", "rb") as f:
        geo_df = gpd.read_file(f)

    transformations = get_transformation(sdata[points_key], get_all=True).copy()
    geo_df = ShapesModel.parse(geo_df, transformations=transformations)

    return adata, geo_df
