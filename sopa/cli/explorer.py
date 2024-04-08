from __future__ import annotations

import typer

from .utils import SDATA_HELPER

app_explorer = typer.Typer()

PIXELSIZE_DEPRECATED = (
    "`pixelsize` is deprecated and will be removed in future versions. Use `pixel_size` instead."
)


@app_explorer.command()
def write(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    output_path: str = typer.Option(
        None,
        help="Path to a directory where Xenium Explorer's outputs will be saved. By default, writes to the same path as `sdata_path` but with the `.explorer` suffix",
    ),
    gene_column: str = typer.Option(
        None, help="Column name of the points dataframe containing the gene names"
    ),
    shapes_key: str = typer.Option(
        None,
        help="Sdata key for the boundaries. By default, uses the baysor boundaires, else the cellpose boundaries",
    ),
    pixel_size: float = typer.Option(
        0.2125,
        help="Number of microns in a pixel. Invalid value can lead to inconsistent scales in the Explorer.",
    ),
    pixelsize: float = typer.Option(
        None,
        help=PIXELSIZE_DEPRECATED,
    ),
    lazy: bool = typer.Option(
        True,
        help="If `True`, will not load the full images in memory (except if the image memory is below `ram_threshold_gb`)",
    ),
    ram_threshold_gb: int = typer.Option(
        4,
        help="Threshold (in gygabytes) from which image can be loaded in memory. If `None`, the image is never loaded in memory",
    ),
    mode: str = typer.Option(
        None,
        help="String that indicated which files should be created. `'-ib'` means everything except images and boundaries, while `'+tocm'` means only transcripts/observations/counts/metadata (each letter corresponds to one explorer file). By default, keeps everything",
    ),
    save_h5ad: bool = typer.Option(
        True,
        help="Whether to save the adata as h5ad in the explorer directory (for convenience only, since h5ad is faster to open than the original .zarr table)",
    ),
):
    """Convert a spatialdata object to Xenium Explorer's inputs"""
    import logging
    from pathlib import Path

    from sopa.io.explorer import write
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    if output_path is None:
        output_path = Path(sdata_path).with_suffix(".explorer")

    if pixelsize is not None:
        log = logging.getLogger(__name__)
        log.critical(PIXELSIZE_DEPRECATED)
        pixel_size = pixelsize

    write(
        output_path,
        sdata,
        shapes_key=shapes_key,
        gene_column=gene_column,
        pixel_size=pixel_size,
        lazy=lazy,
        ram_threshold_gb=ram_threshold_gb,
        mode=mode,
        save_h5ad=save_h5ad,
    )


@app_explorer.command()
def update_obs(
    adata_path: str = typer.Argument(
        help="Path to the anndata file (`zarr` or `h5ad`) containing the new observations"
    ),
    output_path: str = typer.Argument(
        help="Path to the Xenium Explorer directory (it will update `analysis.zarr.zip`)",
    ),
):
    """Update the cell categories for the Xenium Explorer's (i.e. what's in `adata.obs`). This is useful when you perform analysis and update your `AnnData` object

    Usage:
        Make sure you have already run `sopa explorer write` before. This command should only be used if you updated `adata.obs`
    """
    from pathlib import Path

    import anndata

    from sopa.io.explorer import write_cell_categories

    path = Path(adata_path)

    if path.is_dir():
        adata = anndata.read_zarr(path)
    else:
        adata = anndata.read_h5ad(path)

    write_cell_categories(output_path, adata)


@app_explorer.command()
def add_aligned(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    image_path: str = typer.Argument(
        help="Path to the image file to be added (`.ome.tif` used in the explorer during alignment)"
    ),
    transformation_matrix_path: str = typer.Argument(
        help="Path to the `matrix.csv` file returned by the Explorer after alignment"
    ),
    original_image_key: str = typer.Option(
        None,
        help="Optional original-image key (of sdata.images) on which the new image will be aligned. This doesn't need to be provided if there is only one image",
    ),
    overwrite: bool = typer.Option(False, help="Whether to overwrite the image if existing"),
):
    """After alignment on the Xenium Explorer, add an image to the SpatialData object"""
    import spatialdata

    from sopa import io
    from sopa.io.explorer.images import align

    sdata = spatialdata.read_zarr(sdata_path)
    image = io.ome_tif(image_path, as_image=True)

    align(
        sdata, image, transformation_matrix_path, overwrite=overwrite, image_key=original_image_key
    )
