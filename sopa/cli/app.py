import ast

import typer

from .annotate import app_annotate
from .check import app_check
from .patchify import app_patchify
from .resolve import app_resolve
from .segmentation import app_segmentation

option = typer.Option()

app = typer.Typer()
app.add_typer(
    app_annotate,
    name="annotate",
    help="Perform cell-type annotation (based on transcripts and/or channel intensities)",
)
app.add_typer(
    app_segmentation,
    name="segmentation",
    help="Perform cell segmentation on patches (you first need to run 'sopa patchify'). NB: for 'baysor', use directly the 'baysor' command line.",
)
app.add_typer(
    app_resolve, name="resolve", help="Resolve the segmentation conflicts over patches overlaps"
)
app.add_typer(
    app_patchify,
    name="patchify",
    help="Create patches with overlaps. Afterwards, segmentation will be run on each patch",
)
app.add_typer(
    app_check,
    name="check",
    help="Run some sanity checks (e.g., on the YAML config, on the tangram reference, ...)",
)


@app.command()
def read(
    data_path: str,
    technology: str = None,
    sdata_path: str = None,
    config_path: str = None,
    kwargs: str = typer.Option(default={}, callback=ast.literal_eval),
):
    """Read any technology data, and write a standardized SpatialData object

    Args:\n
        data_path: Path to one data sample (most of the time, this corresponds to a directory)\n
        technology: Name of the technology used to collected the data (e.g., 'xenium', 'merfish', ...)\n
        sdata_path: Optional path to write the SpatialData object. If not provided, will write to the '{data_path}.zarr' directory\n
        config_path: Path to the snakemake config. This can be useful in order not to provide the 'technology' and the 'kwargs' arguments\n
        kwargs: Dictionary provided to the reader function.\n
    """
    from pathlib import Path

    from sopa import io

    sdata_path = Path(data_path).with_suffix(".zarr") if sdata_path is None else sdata_path

    assert (
        technology is not None or config_path is not None
    ), "Provide the argument `--technology` or `--config-path`"

    if config_path is not None:
        assert not kwargs, "Provide either a path to a config, or some kwargs, but not both"
        with open(config_path, "r") as f:
            import yaml

            config = yaml.safe_load(f)

        technology = config["read"]["technology"]
        kwargs = config["read"]["kwargs"]

    assert hasattr(
        io, technology
    ), f"Technology {technology} unknown. Currently available: xenium, merscope, cosmx, qptiff"

    sdata = getattr(io, technology)(data_path, **kwargs)
    io.write_standardized(sdata, sdata_path, delete_table=True)


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
    and transfer back the resulting polygon

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

    from sopa.io.standardize import read_zarr_standardized
    from sopa.utils.polygon_crop import intermediate_selection, polygon_selection

    if sdata_path is None:
        assert (
            intermediate_image is not None and intermediate_polygon is not None
        ), "When no --sdata_path is provided, both --intermediate_image and --intermediate_polygon have to be provided"

        intermediate_selection(intermediate_image, intermediate_polygon, margin_ratio)
        return

    sdata = read_zarr_standardized(sdata_path)

    polygon_selection(
        sdata, intermediate_image, intermediate_polygon, list(channels), scale_factor, margin_ratio
    )


@app.command()
def aggregate(
    sdata_path: str,
    gene_column: str = None,
    average_intensities: bool = False,
    min_transcripts: int = 0,
    min_intensity_ratio: float = 0,
):
    """Create an `anndata` table containing the transcript count and/or the channel intensities per cell

    Args:\n
        sdata_path: Path to the SpatialData zarr directory\n
        gene_column: Column of the transcript dataframe representing the gene names. If not provided, it will not compute transcript count\n
        average_intensities: Whether to average the channel intensities inside each cell\n
        min_transcripts: Cells with less transcript than this integer will be filtered\n
        min_intensity_ratio: Cells whose mean channel intensity is less than min_intensity_ratio * quantile_90 will be filtered\n
    """
    from sopa.io.standardize import read_zarr_standardized
    from sopa.segmentation.aggregate import Aggregator

    sdata = read_zarr_standardized(sdata_path)

    aggregator = Aggregator(sdata)
    aggregator.update_table(gene_column, average_intensities, min_transcripts, min_intensity_ratio)


@app.command()
def report(
    sdata_path: str,
    path: str,
):
    """Create a HTML report of the pipeline run and some quality controls

    Args:\n
        sdata_path: Path to the SpatialData zarr directory\n
        path: Path to the HTML report\n
    """
    from sopa.io.report import write_report
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    write_report(path, sdata)


@app.command()
def explorer(
    sdata_path: str,
    path: str,
    gene_column: str = None,
    shapes_key: str = None,
    lazy: bool = True,
    ram_threshold_gb: int = None,
    save_image_mode: int = 1,
):
    """Convert a spatialdata object to Xenium Explorer's inputs

    Args:\n
        sdata_path: Path to the SpatialData zarr directory\n
        path: Path to a directory where Xenium Explorer's outputs will be saved\n
        shapes_key: Key for the boundaries. By default, uses the baysor boundaires, else the cellpose boundaries.\n
        gene_column: Column name of the points dataframe containing the gene names.\n
        lazy: If `True`, will not load the full images in memory (except if the image memory is below `ram_threshold_gb`).\n
        ram_threshold_gb: Threshold (in gygabytes) from which image can be loaded in memory. If `None`, the image is never loaded in memory.\n
        save_image_mode: `1` is normal mode. `0` doesn't save the image. `2` saves **only** the image.\n
    """
    from sopa.io.explorer import write_explorer
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)
    write_explorer(
        path,
        sdata,
        shapes_key=shapes_key,
        gene_column=gene_column,
        lazy=lazy,
        ram_threshold_gb=ram_threshold_gb,
        save_image_mode=save_image_mode,
    )
