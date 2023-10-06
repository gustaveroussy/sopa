import ast

import typer

from .annotate import app_annotate
from .patchify import app_patchify
from .resolve import app_resolve
from .segmentation import app_segmentation

option = typer.Option()

app = typer.Typer()
app.add_typer(app_annotate, name="annotate")
app.add_typer(app_segmentation, name="segmentation")
app.add_typer(app_resolve, name="resolve")
app.add_typer(app_patchify, name="patchify")


@app.command()
def read(
    data_path: str,
    technology: str = option,
    sdata_path: str = None,
    config_path: str = None,
    kwargs: str = typer.Option(default={}, callback=ast.literal_eval),
):
    from pathlib import Path

    from sopa import io

    sdata_path = Path(data_path).with_suffix(".zarr") if sdata_path is None else sdata_path

    assert hasattr(
        io, technology
    ), f"Technology {technology} unknown. Currently available: xenium, merscope, cosmx, qptiff"

    if config_path is not None:
        assert not kwargs, "Provide either a path to a config, or some kwargs, but not both"
        with open("data.yaml", "r") as f:
            import yaml

            kwargs = yaml.safe_load(f)["reader"]["kwargs"]

    sdata = getattr(io, technology)(data_path, **kwargs)
    io.write_standardized(sdata, sdata_path)


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
def aggregate(
    sdata_path: str,
    intensity_mean: bool = True,
    gene_column: str = None,
):
    import spatialdata

    from sopa.segmentation.aggregate import aggregate

    sdata = spatialdata.read_zarr(sdata_path)

    aggregate(sdata, gene_column, intensity_mean)


@app.command()
def explorer(sdata_path: str, path: str, shapes_key: str = None, gene_column: str = None):
    """Convert a spatialdata object to Xenium Explorer's inputs

    Args:\n
        sdata_path: Path to the sdata.zarr directory\n
        path: Path to a directory where Xenium Explorer's outputs will be saved\n
        shapes_key: Key for the boundaries\n
    """
    from spatialdata import SpatialData

    from sopa.io.explorer import write_explorer

    sdata = SpatialData.read(sdata_path)
    write_explorer(
        path, sdata, shapes_key=shapes_key, gene_column=gene_column
    )  # TODO: add more args
