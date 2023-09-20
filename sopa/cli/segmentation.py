import typer

app_segmentation = typer.Typer()
option = typer.Option()


@app_segmentation.command()
def cellpose(
    sdata_path: str,
    diameter: float = option,
    channels: list[str] = option,
    tile_width: int = option,
    tile_overlap: int = option,
    expand_radius: float = option,
):
    """Perform cellpose segmentation on a SpatialData object

    Args:\n
        sdata_path: Path to the sdata.zarr directory.\n
        name: Name of the segmentation method. Can be either 'cellpose' or 'baysor'.\n
    """
    import spatialdata

    from sopa.segmentation import StainingSegmentation, shapes
    from sopa.segmentation.cellpose import cellpose_patch
    from sopa.segmentation.update import update

    sdata = spatialdata.read_zarr(sdata_path)

    method = cellpose_patch(diameter, channels)
    segmentation = StainingSegmentation(sdata, method, channels, tile_width, tile_overlap)
    polygons = segmentation.run_patches()
    polygons = shapes.expand(polygons, expand_radius)

    update(sdata, polygons, segmentation.image_key)


@app_segmentation.command("baysor")
def test():
    ...
