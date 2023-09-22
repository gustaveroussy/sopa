import typer

app_segmentation = typer.Typer()
option = typer.Option()


@app_segmentation.command()
def cellpose(
    sdata_path: str,
    diameter: float = option,
    channels: list[str] = option,
    tile_width: int = typer.Option(default=None),
    tile_overlap: int = typer.Option(default=None),
    expand_radius: int = typer.Option(default=0),
    patch_index: int = typer.Option(default=None),
    patch_dir: str = typer.Option(default=None),
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
    segmentation = StainingSegmentation(sdata, method, channels)

    if patch_index is not None:
        segmentation.write_patch_polygons(patch_dir, patch_index)
        return

    polygons = segmentation.run_patches(tile_width, tile_overlap)
    polygons = shapes.expand(polygons, expand_radius)

    update(sdata, polygons, segmentation.image_key)


@app_segmentation.command()
def resolve(
    sdata_path: str,
    patch_dir: str = option,
    expand_radius: float = typer.Option(default=0),
):
    from pathlib import Path

    import spatialdata
    import zarr
    from shapely.geometry import Polygon
    from tqdm import tqdm

    from sopa.segmentation import shapes
    from sopa.segmentation.update import update
    from sopa.utils.utils import _get_spatial_image

    sdata = spatialdata.read_zarr(sdata_path)

    image_key, _ = _get_spatial_image(sdata)

    polygons = []

    files = [f for f in Path(patch_dir).iterdir() if f.suffix == ".zip"]
    for file in tqdm(files):
        z = zarr.open(file, mode="r")
        for _, coords_zarr in z.arrays():
            polygons.append(Polygon(coords_zarr[:]))
    polygons = shapes.solve_conflicts(polygons)
    polygons = shapes.expand(polygons, expand_radius)

    update(sdata, polygons, image_key)
