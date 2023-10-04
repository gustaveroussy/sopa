import typer

app_segmentation = typer.Typer()
option = typer.Option()


@app_segmentation.command()
def cellpose(
    sdata_path: str,
    diameter: float = option,
    channels: list[str] = option,
    flow_threshold: int = option,
    cellprob_threshold: int = option,
    patch_width: int = typer.Option(default=None),
    patch_overlap: int = typer.Option(default=None),
    expand_radius: int = typer.Option(default=0),
    patch_index: int = typer.Option(default=None),
    patch_dir: str = typer.Option(default=None),
):
    """Perform cellpose segmentation on a SpatialData object

    Args:\n
        sdata_path: Path to the sdata.zarr directory.\n
    """
    import spatialdata

    from sopa.segmentation import shapes
    from sopa.segmentation.cellpose import cellpose_patch
    from sopa.segmentation.cellpose.update import add_shapes
    from sopa.segmentation.stainings import StainingSegmentation

    sdata = spatialdata.read_zarr(sdata_path)

    method = cellpose_patch(
        diameter, channels, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold
    )
    segmentation = StainingSegmentation(sdata, method, channels)

    if patch_index is not None:
        segmentation.write_patch_cells(patch_dir, patch_index)
        return

    cells = segmentation.run_patches(patch_width, patch_overlap)
    cells = shapes.expand(cells, expand_radius)

    add_shapes(sdata, cells, segmentation.image_key)
