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
    tile_width: int = typer.Option(default=None),
    tile_overlap: int = typer.Option(default=None),
    expand_radius: int = typer.Option(default=0),
    patch_index: int = typer.Option(default=None),
    patch_dir: str = typer.Option(default=None),
):
    """Perform cellpose segmentation on a SpatialData object

    Args:\n
        sdata_path: Path to the sdata.zarr directory.\n
    """
    import spatialdata

    from sopa.segmentation import StainingSegmentation, shapes
    from sopa.segmentation.cellpose import cellpose_patch
    from sopa.segmentation.update import update

    sdata = spatialdata.read_zarr(sdata_path)

    method = cellpose_patch(
        diameter, channels, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold
    )
    segmentation = StainingSegmentation(sdata, method, channels)

    if patch_index is not None:
        segmentation.write_patch_polygons(patch_dir, patch_index)
        return

    polygons = segmentation.run_patches(tile_width, tile_overlap)
    polygons = shapes.expand(polygons, expand_radius)

    update(sdata, polygons, segmentation.image_key)
