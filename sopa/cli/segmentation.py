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
    model_type: str = "cyto2",
    patch_width: int = typer.Option(default=None),
    patch_overlap: int = typer.Option(default=None),
    expand_radius: int = typer.Option(default=0),
    patch_index: int = typer.Option(default=None),
    patch_dir: str = typer.Option(default=None),
):
    """Perform cellpose segmentation. This can be done on all patches directly, or on one individual patch (provide 'patch_dir' and 'patch_index')

    Args:\n
        sdata_path: Path to the SpatialData zarr directory\n
        diameter: Cellpose diameter parameter\n
        channels: Names of the channels used for Cellpose. If one channel, then provide just a nucleus channel. If two channels, this is the nucleus and then the cytoplasm channel.\n
        flow_threshold: Cellpose flow_threshold parameter\n
        cellprob_threshold: Cellpose cellprob_threshold parameter\n
        model_type: Name of the cellpose model\n
        patch_width: Ignore this if you already run 'sopa patchify'. Patch width in pixels.\n
        patch_overlap: Ignore this if you already run 'sopa patchify'. Patches overlaps in pixels.\n
        expand_radius: Ignore this if you already run 'sopa patchify'. Cell boundaries radius expansion in pixels.\n
        patch_index: Index of the patch on which cellpose should be run. NB: the number of patches is `len(sdata['sopa_patches'])`.\n
        patch_dir: Path to the temporary cellpose directory inside which we will store each individual patch segmentation\n
    """
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation import shapes
    from sopa.segmentation.cellpose import cellpose_patch
    from sopa.segmentation.cellpose.update import add_shapes
    from sopa.segmentation.stainings import StainingSegmentation

    sdata = read_zarr_standardized(sdata_path)

    method = cellpose_patch(
        diameter,
        channels,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        model_type=model_type,
    )
    segmentation = StainingSegmentation(sdata, method, channels)

    if patch_index is not None:
        segmentation.write_patch_cells(patch_dir, patch_index)
        return

    cells = segmentation.run_patches(patch_width, patch_overlap)
    cells = shapes.expand(cells, expand_radius)

    add_shapes(sdata, cells, segmentation.image_key)
