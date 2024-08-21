from __future__ import annotations

import ast

import typer

from .._constants import SopaKeys
from .utils import SDATA_HELPER

app_patchify = typer.Typer()


@app_patchify.command()
def image(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_width_pixel: float = typer.Option(5000, help="Width (and height) of each patch in pixels"),
    patch_overlap_pixel: float = typer.Option(
        100,
        help="Number of overlapping pixels between the patches. We advise to choose approximately twice the diameter of a cell",
    ),
):
    """Prepare patches for staining-based segmentation (including Cellpose)"""
    from sopa._constants import SopaFiles
    from sopa._sdata import get_spatial_image
    from sopa.io.standardize import read_zarr_standardized
    from sopa.patches import Patches2D

    sdata = read_zarr_standardized(sdata_path)

    image_key, _ = get_spatial_image(sdata, return_key=True)

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
        help="Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation",
    ),
    unassigned_value: str = typer.Option(
        None,
        help="If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)",
    ),
    use_prior: bool = typer.Option(
        False,
        help="Whether to use cellpose segmentation as a prior for baysor (if True, make sure to first run Cellpose)",
    ),
):
    """Prepare patches for transcript-based segmentation with Baysor"""
    from sopa._constants import SopaFiles, SopaKeys

    from .utils import _default_boundary_dir

    if baysor_temp_dir is None:
        baysor_temp_dir = _default_boundary_dir(sdata_path, SopaKeys.BAYSOR_BOUNDARIES)
    return _patchify_transcripts(
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
        help="Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation.",
    ),
    unassigned_value: str = typer.Option(
        None,
        help="If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)",
    ),
    min_transcripts_per_patch: int = typer.Option(1000, help="Minimum number of transcripts per patch"),
    min_cells_per_patch: int = typer.Option(1, help="Minimum number of cells per patch"),
    shapes_key: str = typer.Option(
        SopaKeys.CELLPOSE_BOUNDARIES,
        help="Name of the nuclei boundaries shape use for the prior and centroid in the sdata object",
    ),
    use_prior: bool = typer.Option(
        False,
        help="Whether to use cellpose segmentation as a prior for baysor (if True, make sure to first run Cellpose)",
    ),
):
    """Prepare patches for transcript-based segmentation with ComSeg"""

    from sopa._constants import SopaFiles, SopaKeys

    from .utils import _default_boundary_dir

    if comseg_temp_dir is None:
        comseg_temp_dir = _default_boundary_dir(sdata_path, SopaKeys.COMSEG_BOUNDARIES)

    return _patchify_transcripts(
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
        use_prior=use_prior or cell_key is None,
        min_transcripts_per_patch=min_transcripts_per_patch,
        min_cells_per_patch=min_cells_per_patch,
        shapes_key=shapes_key,
    )


def _patchify_transcripts(
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
    unassigned_value: str,
    min_transcripts_per_patch: int = 4000,
    min_cells_per_patch: int | None = None,
    shapes_key: str = SopaKeys.CELLPOSE_BOUNDARIES,
):
    """Prepare patches for transcript-based segmentation for the different available methods (Baysor, ComSeg).

    Args:
        sdata_path: Path to the SpatialData object
        patch_width_microns: Width (and height) of each patch in microns
        patch_overlap_microns: Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell
        temp_dir: Temporary directory where baysor inputs and outputs will be saved. By default, uses `.sopa_cache/baysor_boundaries`
        filename: Name of the file to indicating the patch's index
        config_name: Name of the config file created for each patch
        config_path: Path to the config file (you can also directly provide the argument via the `config` option)
        config: Dictionnary of parameters
        cell_key: Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation.
        use_prior: Whether to use cellpose segmentation as a prior for baysor and comseg (if True, make sure to first run Cellpose)
        unassigned_value: If cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)

    """
    from sopa._constants import SopaFiles
    from sopa._sdata import get_spatial_element
    from sopa.io.standardize import read_zarr_standardized
    from sopa.patches import Patches2D

    sdata = read_zarr_standardized(sdata_path)

    if isinstance(unassigned_value, str) and unassigned_value.isdigit():
        unassigned_value = int(unassigned_value)

    assert (
        config or config_path is not None
    ), "Provide '--config-path', the path to a Baysor config file (toml) or comseg file (jsons)"

    df_key, _ = get_spatial_element(sdata.points, return_key=True)
    patches = Patches2D(sdata, df_key, patch_width_microns, patch_overlap_microns)
    valid_indices = patches.patchify_transcripts(
        temp_dir,
        cell_key,
        unassigned_value,
        use_prior,
        config,
        config_path,
        config_name=config_name,
        min_transcripts_per_patch=min_transcripts_per_patch,
        shapes_key=shapes_key,
    )

    if filename == SopaFiles.PATCHES_DIRS_COMSEG:
        assert min_cells_per_patch is not None, "For ComSeg, you must provide min_cells_per_patch"

        valid_indices_centoid = patches.patchify_centroids(
            temp_dir,
            shapes_key=shapes_key,
            min_cells_per_patch=min_cells_per_patch,
            cell_key=cell_key,
        )
        valid_indices = list(set(valid_indices_centoid).intersection(set(valid_indices)))

    _save_cache(sdata_path, filename, "\n".join(map(str, valid_indices)))


def _save_cache(sdata_path: str, filename: str, content: str):
    from pathlib import Path

    from sopa._constants import SopaFiles

    cache_file = Path(sdata_path) / SopaFiles.SOPA_CACHE_DIR / filename
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    with open(cache_file, "w") as f:
        f.write(content)
