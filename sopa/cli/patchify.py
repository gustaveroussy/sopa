from __future__ import annotations

import ast

import typer

from .utils import SDATA_HELPER

app_patchify = typer.Typer()


@app_patchify.command()
def image(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_width_pixel: float = typer.Option(
        5000, help="Width (and height) of each patch in pixels"
    ),
    patch_overlap_pixel: float = typer.Option(
        100,
        help="Number of overlapping pixels between the patches. We advise to choose approximately twice the diameter of a cell",
    ),
):
    """Prepare patches for staining-based segmentation (including Cellpose)"""
    from sopa._constants import SopaFiles
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized, sanity_check
    from sopa.segmentation.patching import Patches2D

    sdata = read_zarr_standardized(sdata_path)
    sanity_check(sdata)

    image_key = get_key(sdata, "images")

    patches = Patches2D(sdata, image_key, patch_width_pixel, patch_overlap_pixel)
    patches.write()

    _save_cache(sdata_path, SopaFiles.PATCHES_FILE_IMAGE, str(len(patches)))


@app_patchify.command()
def baysor(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_width_microns: float = typer.Option(help="Width (and height) of each patch in microns"),
    patch_overlap_microns: float = typer.Option(
        help="Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell"
    ),
    baysor_temp_dir: str = typer.Option(
        None,
        help="Temporary directory where baysor inputs and outputs will be saved. By default, uses `.sopa_cache/baysor_boundaries`",
    ),
    config_path: str = typer.Option(
        None,
        help="Path to the baysor config (you can also directly provide the argument via the `config` option)",
    ),
    config: str = typer.Option(
        default={},
        callback=ast.literal_eval,
        help="Dictionnary of baysor parameters",
    ),
    cell_key: str = typer.Option(
        None,
        help="Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation",
    ),
    unassigned_value: int = typer.Option(
        None,
        help="If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)",
    ),
    use_prior: bool = typer.Option(
        False,
        help="Whether to use cellpose segmentation as a prior for baysor (if True, make sure to first run Cellpose)",
    ),
):
    """Prepare the patches for Baysor segmentation"""
    from sopa._constants import SopaFiles, SopaKeys
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized, sanity_check
    from sopa.segmentation.patching import Patches2D

    from .utils import _default_boundary_dir

    sdata = read_zarr_standardized(sdata_path)
    sanity_check(sdata)

    assert (
        config or config_path is not None
    ), "Provide '--config-path', the path to a Baysor config file (toml)"

    if baysor_temp_dir is None:
        baysor_temp_dir = _default_boundary_dir(sdata_path, SopaKeys.BAYSOR_BOUNDARIES)

    df_key = get_key(sdata, "points")
    patches = Patches2D(sdata, df_key, patch_width_microns, patch_overlap_microns)
    valid_indices = patches.patchify_transcripts(
        baysor_temp_dir, cell_key, unassigned_value, use_prior, config, config_path
    )

    _save_cache(sdata_path, SopaFiles.PATCHES_DIRS_BAYSOR, "\n".join(map(str, valid_indices)))


def _save_cache(sdata_path: str, filename: str, content: str):
    from pathlib import Path

    from sopa._constants import SopaFiles

    cache_file = Path(sdata_path) / SopaFiles.SOPA_CACHE_DIR / filename
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    with open(cache_file, "w") as f:
        f.write(content)
