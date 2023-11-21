import typer

app_explorer = typer.Typer()
option = typer.Option()


@app_explorer.command()
def write(
    sdata_path: str,
    output_path: str = None,
    gene_column: str = None,
    shapes_key: str = None,
    lazy: bool = True,
    ram_threshold_gb: int = 4,
    mode: str = None,
):
    """Convert a spatialdata object to Xenium Explorer's inputs

    [Args]\n
        sdata_path: Path to the SpatialData zarr directory\n
    \n
    [Options]\n
        output_path: Path to a directory where Xenium Explorer's outputs will be saved. By default, writes to the same path as `sdata_path` but with the `.explorer` suffix\n
        shapes_key: Key for the boundaries. By default, uses the baysor boundaires, else the cellpose boundaries.\n
        gene_column: Column name of the points dataframe containing the gene names.\n
        lazy: If `True`, will not load the full images in memory (except if the image memory is below `ram_threshold_gb`).\n
        ram_threshold_gb: Threshold (in gygabytes) from which image can be loaded in memory. If `None`, the image is never loaded in memory.\n
        mode: string that indicated which files should be created. "-ib" means everything except images and boundaries, while "+tocm" means only transcripts/observations/counts/metadata (each letter corresponds to one explorer file). By default, keeps everything.\n
    """
    from pathlib import Path

    from sopa.io.explorer import write_explorer
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    if output_path is None:
        output_path = Path(sdata_path).with_suffix(".explorer")

    write_explorer(
        output_path,
        sdata,
        shapes_key=shapes_key,
        gene_column=gene_column,
        lazy=lazy,
        ram_threshold_gb=ram_threshold_gb,
        mode=mode,
    )
