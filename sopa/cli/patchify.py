import ast

import typer

app_patchify = typer.Typer()
option = typer.Option()


@app_patchify.command()
def cellpose(
    sdata_path: str,
    patch_width_pixel: float = None,
    patch_overlap_pixel: float = None,
):
    """Prepare patches for Cellpose segmentation

    [Args]\n
        sdata_path: Path to the SpatialData zarr directory\n

    [Options]\n
        patch_width_pixel: Width (and height) of each patch in pixels\n
        patch_overlap_pixel: Number of overlapping pixels between the patches. We advise to choose approximately twice the diameter of a cell\n
    """
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
    sdata_path: str,
    patch_width_microns: float = option,
    patch_overlap_microns: float = option,
    baysor_temp_dir: str = option,
    config: str = typer.Option(default={}, callback=ast.literal_eval),
    cell_key: str = None,
    unassigned_value: int = None,
    use_prior: bool = False,
):
    """Prepare the patches for Baysor segmentation

    [Args]\n
        sdata_path: Path to the SpatialData zarr directory\n
    \n
    [Options]\n
        patch_width_microns: Width (and height) of each patch in microns\n
        patch_overlap_microns: Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell\n
        baysor_temp_dir: Temporary directory where baysor inputs and outputs will be saved\n
        config: Path to the snakemake config containing the baysor arguments\n
        cell_key: Optional column of the transcripts dataframe that indicates in which cell-id each transcript is\n
        unassigned_value: If 'cell_key' is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)\n
        use_prior: Whether to use cellpose segmentation as a prior for baysor\n
    """
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
