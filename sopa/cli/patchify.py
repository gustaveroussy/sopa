import ast

import typer

from .utils import SDATA_HELPER

app_patchify = typer.Typer()


@app_patchify.command()
def cellpose(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_width_pixel: float = typer.Option(
        None, help="Width (and height) of each patch in pixels"
    ),
    patch_overlap_pixel: float = typer.Option(
        None,
        help="Number of overlapping pixels between the patches. We advise to choose approximately twice the diameter of a cell",
    ),
):
    """Prepare patches for Cellpose segmentation"""
    from pathlib import Path

    from sopa._constants import SopaFiles
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized, sanity_check
    from sopa.segmentation.patching import Patches2D

    sdata = read_zarr_standardized(sdata_path)
    sanity_check(sdata)

    image_key = get_key(sdata, "images")

    patches = Patches2D(sdata, image_key, patch_width_pixel, patch_overlap_pixel)
    patches.write()

    with open(Path(sdata_path) / SopaFiles.SMK_DIR / SopaFiles.PATCHES_FILE_CELLPOSE, "w") as f:
        f.write(str(len(patches)))


@app_patchify.command()
def baysor(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_width_microns: float = typer.Option(help="Width (and height) of each patch in microns"),
    patch_overlap_microns: float = typer.Option(
        help="Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell"
    ),
    baysor_temp_dir: str = typer.Option(
        help="Temporary directory where baysor inputs and outputs will be saved"
    ),
    config: str = typer.Option(
        default={},
        callback=ast.literal_eval,
        help="Path to the snakemake config containing the baysor arguments",
    ),
    cell_key: str = typer.Option(
        None,
        help="Optional column of the transcripts dataframe that indicates in which cell-id each transcript is",
    ),
    unassigned_value: int = typer.Option(
        None,
        help="If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)",
    ),
    use_prior: bool = typer.Option(
        False, help="Whether to use cellpose segmentation as a prior for baysor"
    ),
):
    """Prepare the patches for Baysor segmentation"""
    from pathlib import Path

    from sopa._constants import SopaFiles
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized, sanity_check
    from sopa.segmentation.baysor.prepare import to_toml
    from sopa.segmentation.patching import Patches2D

    sdata = read_zarr_standardized(sdata_path)
    sanity_check(sdata)

    assert config is not None

    df_key = get_key(sdata, "points")
    patches = Patches2D(sdata, df_key, patch_width_microns, patch_overlap_microns)
    valid_indices = patches.patchify_transcripts(
        baysor_temp_dir, cell_key, unassigned_value, use_prior
    )

    for i in range(len(patches)):
        path = Path(baysor_temp_dir) / str(i) / SopaFiles.BAYSOR_CONFIG
        to_toml(path, config)

    with open(Path(sdata_path) / SopaFiles.SMK_DIR / SopaFiles.PATCHES_DIRS_BAYSOR, "w") as f:
        f.write("\n".join(map(str, valid_indices)))
