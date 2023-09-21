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
    expand_radius: float = typer.Option(default=0),
    patch_attrs_file: str = typer.Option(default=None),
    patch_file: str = typer.Option(default=None),
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

    assert not (
        patch_attrs_file and patch_file
    ), f"Choose either --patch_file to run one patch, or --patch_attrs_file to prepare the parralel segmentation"

    if patch_attrs_file is not None:
        from sopa.utils.tiling import Tiles2D
        from sopa.utils.utils import _get_spatial_image

        image = _get_spatial_image(sdata)

        tiles = Tiles2D.from_image(image, tile_width, tile_overlap)
        tiles.write(patch_attrs_file)
        print("Saved patches paths. You can now run segmentation on each patch on separate jobs.")
        return

    method = cellpose_patch(diameter, channels)
    segmentation = StainingSegmentation(sdata, method, channels, tile_width, tile_overlap)

    if patch_file is not None:
        segmentation.write_patch_polygons(patch_file)
        return

    polygons = segmentation.run_patches()
    polygons = shapes.expand(polygons, expand_radius)

    update(sdata, polygons, segmentation.image_key)


@app_segmentation.command()
def resolve(
    sdata_path: str,
    patch_attrs_file: str = option,
    expand_radius: float = typer.Option(default=0),
):
    import spatialdata
    import zarr
    from shapely.geometry import Polygon
    from tqdm import tqdm

    from sopa.segmentation import shapes
    from sopa.segmentation.update import update
    from sopa.utils.utils import _get_spatial_image

    sdata = spatialdata.read_zarr(sdata_path)

    image_key, _ = _get_spatial_image(sdata)

    with open(patch_attrs_file, "r") as f:
        patches_paths = [line.rstrip() for line in f]

    polygons = []

    for patches_path in tqdm(patches_paths):
        z = zarr.open(patches_path, mode="r")
        for _, coords_zarr in z.arrays():
            polygons.append(Polygon(coords_zarr[:]))

    polygons = shapes.solve_conflicts(polygons)
    polygons = shapes.expand(polygons, expand_radius)

    update(sdata, polygons, image_key)
