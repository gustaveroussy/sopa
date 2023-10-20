import typer

app_resolve = typer.Typer()
option = typer.Option()


@app_resolve.command()
def cellpose(
    sdata_path: str,
    patch_dir: str = option,
    expand_radius: float = typer.Option(default=0),
):
    import spatialdata

    from sopa._sdata import get_key
    from sopa.segmentation import shapes
    from sopa.segmentation.cellpose.update import add_shapes
    from sopa.segmentation.stainings import StainingSegmentation

    sdata = spatialdata.read_zarr(sdata_path)

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
    import spatialdata

    from sopa.segmentation.baysor.resolve import resolve

    assert (
        baysor_temp_dir is not None or patches_dirs is not None
    ), "Provide either a baysor directory (--baysor_temp_dir) or a list of all subdirectories (--patches_dirs)"

    sdata = spatialdata.read_zarr(sdata_path)

    resolve(sdata, baysor_temp_dir, gene_column, patches_dirs, min_area, expand_radius)
