from __future__ import annotations

import typer

from .utils import SDATA_HELPER

app_resolve = typer.Typer()


@app_resolve.command()
def cellpose(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    patch_dir: list[str] = typer.Option(
        [],
        help="Directory containing the cellpose segmentation on patches (or multiple directories if using multi-step segmentation). By default, uses the `.sopa_cache/cellpose_boundaries` directory",
    ),
):
    """Resolve patches conflicts after cellpose segmentation"""
    from sopa._constants import SopaKeys

    from .utils import _default_boundary_dir

    if not len(patch_dir):
        patch_dir = [_default_boundary_dir(sdata_path, SopaKeys.CELLPOSE_BOUNDARIES)]

    _resolve_generic(sdata_path, patch_dir, SopaKeys.CELLPOSE_BOUNDARIES)


@app_resolve.command()
def generic(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    method_name: str = typer.Option(
        help="Name of the method used during segmentation. This is also the key correspnding to the boundaries in `sdata.shapes`"
    ),
    patch_dir: list[str] = typer.Option(
        [],
        help="Directory containing the generic segmentation on patches (or multiple directories if using multi-step segmentation). By default, uses the `.sopa_cache/<method_name>` directory",
    ),
):
    """Resolve patches conflicts after generic segmentation"""
    from .utils import _default_boundary_dir

    if not len(patch_dir):
        patch_dir = [_default_boundary_dir(sdata_path, method_name)]

    _resolve_generic(sdata_path, patch_dir, method_name)


def _resolve_generic(sdata_path: str, patch_dirs: list[str], shapes_key: str):
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation import shapes
    from sopa.segmentation.stainings import StainingSegmentation

    sdata = read_zarr_standardized(sdata_path)

    assert len(patch_dirs) > 0, "No patch directory was provided, cannot load cells"

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
        None,
        help="Path to the directory containing all the baysor patches (see `sopa patchify`). By default, uses the `.sopa_cache/baysor_boundaries` directory",
    ),
    min_area: float = typer.Option(
        0, help="Cells with an area less than this value (in microns^2) will be filtered"
    ),
    patches_dirs: list[str] = typer.Option(
        [], help="List of patches directories inside `baysor_temp_dir`"
    ),
):
    """Resolve patches conflicts after baysor segmentation. Provide either `--baysor-temp-dir` or `--patches-dirs`"""
    from sopa._constants import SopaKeys
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation.baysor.resolve import resolve

    from .utils import _default_boundary_dir

    if not len(patches_dirs) and baysor_temp_dir is None:
        baysor_temp_dir = _default_boundary_dir(sdata_path, SopaKeys.BAYSOR_BOUNDARIES)

    sdata = read_zarr_standardized(sdata_path)

    resolve(sdata, baysor_temp_dir, gene_column, patches_dirs, min_area)
