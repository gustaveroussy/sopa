import typer

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
    import sopa
    from sopa._constants import SopaFiles, SopaKeys
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    sopa.make_image_patches(sdata, patch_width_pixel, patch_overlap_pixel)

    n_patches = len(sdata[SopaKeys.PATCHES])
    _save_cache(sdata_path, SopaFiles.PATCHES_FILE_IMAGE, n_patches)


@app_patchify.command()
def transcripts(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_width_microns: float = typer.Option(5000, help="Width (and height) of each patch in microns"),
    patch_overlap_microns: float = typer.Option(
        100,
        help="Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell",
    ),
    unassigned_value: str = typer.Option(
        None,
        help="If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)",
    ),
    prior_shapes_key: str = typer.Option(
        None,
        help="Name of the boundaries element to use as a segmentation prior",
    ),
    write_cells_centroids: bool = typer.Option(
        False,
        help="Whether to also write the centroids of the cells (must be True for ComSeg)",
    ),
):
    """Prepare patches for transcript-based segmentation (including Baysor)"""

    import sopa
    from sopa._constants import SopaFiles, SopaKeys
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    if isinstance(unassigned_value, str) and unassigned_value.isdigit():
        unassigned_value = int(unassigned_value)

    sopa.make_transcript_patches(
        sdata,
        patch_width_microns,
        patch_overlap_microns,
        unassigned_value=unassigned_value,
        prior_shapes_key=prior_shapes_key,
        write_cells_centroids=write_cells_centroids,
    )

    valid_indices = list(sdata[SopaKeys.TRANSCRIPTS_PATCHES].index)
    _save_cache(sdata_path, SopaFiles.PATCHES_FILE_TRANSCRIPTS, "\n".join(map(str, valid_indices)))


def _save_cache(sdata_path: str, filename: str, content: str):
    from pathlib import Path

    from sopa._constants import SopaFiles

    cache_file = Path(sdata_path) / SopaFiles.SOPA_CACHE_DIR / filename
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    with open(cache_file, "w", newline="\n") as f:
        f.write(f"{content}\n")
