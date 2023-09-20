import typer

app = typer.Typer()


@app.command()
def crop(
    sdata_path: str = None,
    intermediate_image: str = None,
    intermediate_polygon: str = None,
    channels: list[str] = None,
    scale_factor: float = 10,
    margin_ratio: float = 0.1,
):
    """Crop an image based on a user-defined polygon (interactive mode).
    If the interactive mode is not available where the data is stored,
    then we can export an intermediate resized image, then make the selection locally,
    and transfer back the resulting polygon.

    Args:\n
        sdata_path: Path to the sdata.zarr directory. Defaults to None.\n
        intermediate_image: Path to the intermediate image, with a .zip extension. Use this only if the interactive mode is not available. Defaults to None.\n
        intermediate_polygon: Path to the intermediate polygon, with a .zip extension. Use this locally, after downloading the intermediate_image. Defaults to None.\n
        channels: List of channel names to be displayed. Defaults to None.\n
        scale_factor: Resize the image by this value (high value for a lower memory usage). Defaults to 10.\n
        margin_ratio: Ratio of the image margin on the display (compared to the image size). Defaults to 0.1.\n
    """
    from .utils import _check_zip

    _check_zip([intermediate_image, intermediate_polygon])

    import spatialdata

    from sopa.utils.polygon_crop import intermediate_selection, polygon_selection

    if sdata_path is None:
        assert (
            intermediate_image is not None and intermediate_polygon is not None
        ), "When no --sdata_path is provided, both --intermediate_image and --intermediate_polygon have to be provided"

        intermediate_selection(intermediate_image, intermediate_polygon, margin_ratio)
        return

    sdata = spatialdata.read_zarr(sdata_path)

    polygon_selection(
        sdata, intermediate_image, intermediate_polygon, list(channels), scale_factor, margin_ratio
    )


@app.command()
def test():
    print("test")
