import ast
import logging

import typer

from .annotate import app_annotate
from .explorer import app_explorer
from .patchify import app_patchify
from .resolve import app_resolve
from .segmentation import app_segmentation
from .utils import SDATA_HELPER

log = logging.getLogger(__name__)

app = typer.Typer()
app.add_typer(
    app_annotate,
    name="annotate",
    help="Perform cell-type annotation (based on transcripts and/or channel intensities)",
)
app.add_typer(
    app_explorer,
    name="explorer",
    help="Convertion to the Xenium Explorer's inputs, and image alignment",
)
app.add_typer(
    app_segmentation,
    name="segmentation",
    help="Perform cell segmentation on patches. NB: for `baysor`, use directly the `baysor` command line.",
)
app.add_typer(app_resolve, name="resolve", help="Resolve the segmentation conflicts over patches overlaps")
app.add_typer(
    app_patchify,
    name="patchify",
    help="Create patches with overlaps. Afterwards, segmentation will be run on each patch",
)


@app.command()
def read(data_path: str = typer.Argument(), technology: str = typer.Option()):
    """Deprecated and will be removed in sopa==2.1.0. Use `sopa convert` instead."""
    raise NameError("`sopa read` is deprecated. Use `sopa convert` instead.")


@app.command()
def convert(
    data_path: str = typer.Argument(
        help="Path to one data sample (most of the time, this corresponds to a directory with images files and eventually a transcript file)"
    ),
    technology: str = typer.Option(
        None,
        help="Name of the technology used to collected the data (`xenium`/`merfish`/`cosmx`/`phenocycler`/`macsima`/`hyperion`)",
    ),
    sdata_path: str = typer.Option(
        None,
        help="Optional path to write the SpatialData object. If not provided, will write to the `{data_path}.zarr` directory",
    ),
    config_path: str = typer.Option(
        None,
        help="Path to the snakemake config. This can be useful in order not to provide the `--technology` and the `--kwargs` arguments",
    ),
    kwargs: str = typer.Option(
        default={},
        callback=ast.literal_eval,
        help="Dictionary provided to the reader function as kwargs",
    ),
):
    """Read any technology data, and write a standardized SpatialData object.

    Either `--technology` or `--config-path` has to be provided."""
    from pathlib import Path

    from sopa import io

    sdata_path = Path(data_path).with_suffix(".zarr") if sdata_path is None else Path(sdata_path)

    io.standardize._check_can_write_zarr(sdata_path)

    assert technology is not None or config_path is not None, "Provide the argument `--technology` or `--config-path`"

    if config_path is not None:
        assert not kwargs, "Provide either a path to a config, or some kwargs, but not both"
        with open(config_path, "r") as f:
            import yaml

            config = yaml.safe_load(f)

        technology = config["read"]["technology"]
        kwargs = config["read"]["kwargs"]

    assert hasattr(
        io, technology
    ), f"Technology {technology} unknown. Currently available: xenium, merscope, visium_hd, cosmx, phenocycler, hyperion, macsima"

    sdata = getattr(io, technology)(data_path, **kwargs)
    io.write_standardized(sdata, sdata_path, delete_table=True)


@app.command()
def crop(
    sdata_path: str = typer.Option(None, help=SDATA_HELPER),
    intermediate_image: str = typer.Option(
        None,
        help="Path to the intermediate image, with a `.zip` extension. Use this only if the interactive mode is not available",
    ),
    intermediate_polygon: str = typer.Option(
        None,
        help="Path to the intermediate polygon, with a `.zip` extension. Use this locally, after downloading the intermediate_image",
    ),
    channels: list[str] = typer.Option(
        None,
        help="List of channel names to be displayed. Optional if there are already only 1 or 3 channels",
    ),
    scale_factor: float = typer.Option(10, help="Resize the image by this value (high value for a lower memory usage)"),
    margin_ratio: float = typer.Option(
        0.1, help="Ratio of the image margin on the display (compared to the image size)"
    ),
):
    """Crop an image based on a user-defined polygon (interactive mode).

    Warning:
        This command is deprecated. Using `napari-spatialdata` instead.
        Provide the name `"region_of_interest"` to your selected ROI.

    Usage:

        - [Locally] Only `--sdata-path` is required

        - [On a cluster] Run `sopa crop` with `--sdata-path` and `--intermediate-image` on the cluster. Then, download the image locally, and run `sopa crop` with `--intermediate-image` and `--intermediate-polygon`. Then, upload the polygon and run `sopa crop` on the cluster with `--sdata-path` and `--intermediate-polygon`.
    """
    from .utils import _check_zip

    _check_zip([intermediate_image, intermediate_polygon])

    from sopa.io.standardize import read_zarr_standardized
    from sopa.utils.crop import intermediate_selection, polygon_selection

    if intermediate_image and intermediate_polygon:
        assert (
            sdata_path is None
        ), "When both --intermediate_image and --intermediate_polygon, sdata_path should not to be provided"

        intermediate_selection(intermediate_image, intermediate_polygon, margin_ratio)
        return

    sdata = read_zarr_standardized(sdata_path)

    polygon_selection(
        sdata,
        intermediate_image,
        intermediate_polygon,
        None if channels is None else list(channels),
        scale_factor,
        margin_ratio,
    )


@app.command()
def aggregate(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    aggregate_genes: bool = typer.Option(None, help="Whether to aggregate the genes (counts) inside each cell"),
    aggregate_channels: bool = typer.Option(
        False, help="Whether to aggregate the channels (intensity) inside each cell"
    ),
    expand_radius_ratio: float = typer.Option(
        default=None,
        help="Cells polygons will be expanded by `expand_radius_ratio * mean_radius` for channels averaging **only**. This help better aggregate boundary stainings",
    ),
    min_transcripts: int = typer.Option(0, help="Cells with less transcript than this integer will be filtered"),
    min_intensity_ratio: float = typer.Option(
        0,
        help="Cells whose mean channel intensity is less than `min_intensity_ratio * quantile_90` will be filtered",
    ),
    image_key: str = typer.Option(
        None,
        help="Optional image key of the SpatialData object. By default, considers the only one image. It can be useful if another image is added later on",
    ),
    method_name: str = typer.Option(
        None, help="If segmentation was performed with a generic method, this is the name of the method used."
    ),
    average_intensities: bool = typer.Option(False, help="[Deprecated] Use `aggregate_channels` instead."),
    gene_column: str = typer.Option(None, help="[Deprecated] Use `aggregate_genes` instead."),
):
    """Create an `anndata` table containing the transcript count and/or the channel intensities per cell"""
    import sopa
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    if gene_column is not None:
        log.warning(
            "The `gene_column` argument is deprecated and will be removed in sopa==2.1.0. Use `aggregate_genes` instead."
        )
        aggregate_genes = True

    if average_intensities:
        log.warning(
            "The `average_intensities` argument is deprecated and will be removed in sopa==2.1.0. Use `aggregate_channels` instead."
        )
        aggregate_channels = True

    sopa.aggregate(
        sdata,
        aggregate_genes=aggregate_genes,
        aggregate_channels=aggregate_channels,
        image_key=image_key,
        shapes_key=method_name,
        min_transcripts=min_transcripts,
        expand_radius_ratio=expand_radius_ratio,
        min_intensity_ratio=min_intensity_ratio,
    )


@app.command()
def report(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    path: str = typer.Argument(help="Path to the HTML report"),
    table_key: str = typer.Option(
        "table", help="Key of the table in the `SpatialData` object to be used for the report"
    ),
):
    """Create a HTML report of the pipeline run and some quality controls"""
    from sopa.io.report import write_report
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    write_report(path, sdata, table_key=table_key)
