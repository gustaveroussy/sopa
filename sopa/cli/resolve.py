import typer

app_resolve = typer.Typer()
option = typer.Option()


@app_resolve.command()
def cellpose(
    sdata_path: str,
    patch_dir: str = option,
    expand_radius: float = typer.Option(default=0),
):
    """Resolve patches conflicts after cellpose segmentation

    [Args]\n
        sdata_path: Path to the SpatialData zarr directory\n
    \n
    [Options]\n
        patch_dir: Directory containing the cellpose segmentation on patches\n
        expand_radius: Number of pixels for radius expansion of each cell boundary\n
    """
    from sopa._sdata import get_key
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation import shapes
    from sopa.segmentation.cellpose.update import add_shapes
    from sopa.segmentation.stainings import StainingSegmentation

    sdata = read_zarr_standardized(sdata_path)

    image_key = get_key(sdata, "images")

    cells = StainingSegmentation.read_patches_cells(patch_dir)
    cells = shapes.solve_conflicts(cells)
    cells = shapes.expand(cells, expand_radius)

    add_shapes(sdata, cells, image_key)


@app_resolve.command()
def baysor(
    sdata_path: str,
    gene_column: str = option,
    baysor_temp_dir: str = None,
    min_area: float = 0,
    expand_radius: float = 0,
    patches_dirs: list[str] = None,
):
    """Resolve patches conflicts after baysor segmentation. Provide either 'baysor_temp_dir' or 'patches_dirs'

    [Args]\n
        sdata_path: Path to the SpatialData zarr directory\n
    \n
    [Options]\n
        gene_column: Column of the transcripts dataframe containing the genes names\n
        baysor_temp_dir: Path to the directory containing all the baysor patches (see 'sopa patchify')\n
        min_area: Cells with an area less than this value (in microns^2) will be filtered\n
        expand_radius: Number of microns for radius expansion of each cell boundary\n
        patches_dirs: List of patches directories inside 'baysor_temp_dir'\n
    """
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation.baysor.resolve import resolve

    assert (
        baysor_temp_dir is not None or patches_dirs is not None
    ), "Provide either a baysor directory (--baysor_temp_dir) or a list of all subdirectories (--patches_dirs)"

    sdata = read_zarr_standardized(sdata_path)

    resolve(sdata, baysor_temp_dir, gene_column, patches_dirs, min_area, expand_radius)
