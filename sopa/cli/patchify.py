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
    from pathlib import Path

    import spatialdata

    sdata = spatialdata.read_zarr(sdata_path)

    from sopa._constants import SopaFiles
    from sopa._sdata import get_key
    from sopa.segmentation.patching import Patches2D

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
    from pathlib import Path

    import spatialdata

    from sopa._constants import SopaFiles
    from sopa._sdata import get_key
    from sopa.segmentation.baysor.prepare import to_toml
    from sopa.segmentation.patching import Patches2D

    sdata = spatialdata.read_zarr(sdata_path)

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
