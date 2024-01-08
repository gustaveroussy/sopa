import typer

from .utils import SDATA_HELPER

app_resolve = typer.Typer()


@app_resolve.command()
def cellpose(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_dir: list[str] = typer.Option(
        help="Directory containing the cellpose segmentation on patches (or multiple directories if using multi-step segmentation)"
    ),
):
    """Resolve patches conflicts after cellpose segmentation"""
    from sopa._constants import SopaKeys

    _resolve_generic(sdata_path, patch_dir, SopaKeys.CELLPOSE_BOUNDARIES)


@app_resolve.command()
def generic(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    method_name: str = typer.Option(
        help="Name of the method used during segmentation. This is also the key correspnding to the boundaries in `sdata.shapes`"
    ),
    patch_dir: list[str] = typer.Option(
        help="Directory containing the generic segmentation on patches (or multiple directories if using multi-step segmentation)"
    ),
):
    """Resolve patches conflicts after generic segmentation"""
    _resolve_generic(sdata_path, patch_dir, method_name)


def _resolve_generic(sdata_path: str, patch_dirs: list[str], shapes_key: str):
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation import shapes
    from sopa.segmentation.stainings import StainingSegmentation

    sdata = read_zarr_standardized(sdata_path)

    image_key = get_key(sdata, "images")

    cells = []
    for patch_dir in patch_dirs:
        cells += StainingSegmentation.read_patches_cells(patch_dir)
    cells = shapes.solve_conflicts(cells)

    StainingSegmentation.add_shapes(sdata, cells, image_key, shapes_key)


@app_resolve.command()
def baysor(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    gene_column: str = typer.Option(
        help="Column of the transcripts dataframe containing the genes names"
    ),
    baysor_temp_dir: str = typer.Option(
        None, help="Path to the directory containing all the baysor patches (see `sopa patchify`)"
    ),
    min_area: float = typer.Option(
        0, help="Cells with an area less than this value (in microns^2) will be filtered"
    ),
    patches_dirs: list[str] = typer.Option(
        None, help="List of patches directories inside `baysor_temp_dir`"
    ),
):
    """Resolve patches conflicts after baysor segmentation. Provide either `--baysor-temp-dir` or `--patches-dirs`"""
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation.baysor.resolve import resolve

    assert (
        baysor_temp_dir is not None or patches_dirs is not None
    ), "Provide either a baysor directory (--baysor_temp_dir) or a list of all subdirectories (--patches_dirs)"

    sdata = read_zarr_standardized(sdata_path)

    resolve(sdata, baysor_temp_dir, gene_column, patches_dirs, min_area)
