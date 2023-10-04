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
    baysor_temp_dir: str = option,
    gene_column: str = option,
    min_area: float = 0,
    n: int = None,
):
    import spatialdata

    from sopa.segmentation.baysor.resolve import resolve

    sdata = spatialdata.read_zarr(sdata_path)

    resolve(sdata, baysor_temp_dir, gene_column, n, min_area)
