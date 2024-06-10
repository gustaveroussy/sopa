from __future__ import annotations

import ast

import typer

from .._constants import SopaKeys
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
    from sopa.patches import Patches2D

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
        help="Dictionnary of baysor parameters, overwrite the config_path argument if provided",
    ),
    cell_key: str = typer.Option(
        None,
        help="Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation"
        f" Default is '{SopaKeys.DEFAULT_CELL_KEY}' if cell_key=None",
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
    """Prepare patches for transcript-based segmentation with baysor"""
    from sopa._constants import SopaKeys, SopaFiles
    from .utils import _default_boundary_dir


    if baysor_temp_dir is None:
        baysor_temp_dir = _default_boundary_dir(sdata_path, SopaKeys.BAYSOR_BOUNDARIES)
    return _transcript_segmentation(
        sdata_path=sdata_path,
        patch_width_microns=patch_width_microns,
        patch_overlap_microns=patch_overlap_microns,
        temp_dir=baysor_temp_dir,
        filename=SopaFiles.PATCHES_DIRS_BAYSOR,
        config_name=SopaFiles.TOML_CONFIG_FILE,
        config_path=config_path,
        config=config,
        cell_key=cell_key,
        unassigned_value=unassigned_value,
        use_prior=use_prior,
    )


@app_patchify.command()
def comseg(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_width_microns: float = typer.Option(help="Width (and height) of each patch in microns"),
    patch_overlap_microns: float = typer.Option(
        help="Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell"
    ),
    comseg_temp_dir: str = typer.Option(
        None,
        help="Temporary directory where baysor inputs and outputs will be saved. By default, uses `.sopa_cache/comseg_boundaries`",
    ),
    config_path: str = typer.Option(
        None,
        help="Path to the ComSeg json config file (you can also directly provide the argument via the `config` option)",
    ),
    config: str = typer.Option(
        default={},
        callback=ast.literal_eval,
        help="Dictionnary of ComSeg parameters, overwrite the config_path argument if provided",
    ),
    cell_key: str = typer.Option(
        None,
        help="Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation."
        f" Default is {SopaKeys.DEFAULT_CELL_KEY} if cell_key=None",
    ),
    unassigned_value: int = typer.Option(
        None,
        help="If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)",
    ),
):
    """Prepare patches for transcript-based segmentation with ComSeg"""
    from sopa._constants import SopaKeys, SopaFiles
    from .utils import _default_boundary_dir

    if comseg_temp_dir is None:
        comseg_temp_dir = _default_boundary_dir(sdata_path, SopaKeys.COMSEG_BOUNDARIES)

    return _transcript_segmentation(
        sdata_path=sdata_path,
        patch_width_microns=patch_width_microns,
        patch_overlap_microns=patch_overlap_microns,
        temp_dir=comseg_temp_dir,
        filename=SopaFiles.PATCHES_DIRS_COMSEG,
        config_name=SopaFiles.JSON_CONFIG_FILE,
        config_path=config_path,
        config=config,
        cell_key=cell_key,
        unassigned_value=unassigned_value,
        use_prior=True,
    )


def _transcript_segmentation(
    sdata_path: str,
    patch_width_microns: float,
    patch_overlap_microns: float,
    temp_dir: str,
    filename: str,
    config_name: str,
    config_path: str,
    config: str,
    cell_key: str,
    use_prior: bool,
    unassigned_value: int,
):
    """Prepare patches for transcript-based segmentation for the different available methods (baysor, comseg)
    args:
        sdata_path (str) : Path to the SpatialData object
        patch_width_microns (float) : Width (and height) of each patch in microns
        patch_overlap_microns (str) : Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell
        temp_dir (str) : Temporary directory where baysor inputs and outputs will be saved. By default, uses `.sopa_cache/baysor_boundaries`"
        filename (str) : Name of the file to indicating the patch's index
        config_name (str) : Name of the config file created for each patch
        config_path (str): "Path to the config file (you can also directly provide the argument via the `config` option)"
        config (str): "Dictionnary of parameters"
        cell_key (str): "Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation.
        " Default is cell if cell_key=None"
        use_prior (bool): "Whether to use cellpose segmentation as a prior for baysor and comseg (if True, make sure to first run Cellpose)"
        unassigned_value (int): "If cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)"

    """
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized, sanity_check
    from sopa.patches import Patches2D

    from .utils import _default_boundary_dir

    sdata = read_zarr_standardized(sdata_path)
    sanity_check(sdata)

    assert (
        config or config_path is not None
    ), "Provide '--config-path', the path to a Baysor config file (toml) or comseg file (jsons)"

    df_key = get_key(sdata, "points")
    patches = Patches2D(sdata, df_key, patch_width_microns, patch_overlap_microns)
    if filename==SopaKeys.COMSEG_BOUNDARIES:
        patches.patchify_centroids(temp_dir)
        assert (
            use_prior
        ), "For ComSeg, you must use the prior segmentation of nuclei or from other staining"
    valid_indices = patches.patchify_transcripts(
        temp_dir,
        cell_key,
        unassigned_value,
        use_prior,
        config,
        config_path,
        config_name=config_name,
    )
    _save_cache(sdata_path, filename, "\n".join(map(str, valid_indices)))


def _save_cache(sdata_path: str, filename: str, content: str):
    from pathlib import Path

    from sopa._constants import SopaFiles

    cache_file = Path(sdata_path) / SopaFiles.SOPA_CACHE_DIR / filename
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    with open(cache_file, "w") as f:
        f.write(content)
