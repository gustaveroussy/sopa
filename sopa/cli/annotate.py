import ast

import typer

app_annotate = typer.Typer()
option = typer.Option()


@app_annotate.command()
def fluorescence(
    sdata_path: str,
    marker_cell_dict: str = typer.Option(default={}, callback=ast.literal_eval),
    key: str = "cell_type",
):
    from pathlib import Path

    from sopa.annotation.fluorescence import higher_z_score
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    assert sdata.table is not None, f"Annotation requires `sdata.table` to be not None"

    higher_z_score(sdata.table, marker_cell_dict, key)
    sdata.table.write_zarr(Path(sdata_path) / "table" / "table")


@app_annotate.command()
def tangram(
    sdata_path: str,
    sc_reference_path: str = option,
    cell_type_key: str = "cell_type",
    bag_size: int = 10_000,
):
    from pathlib import Path

    import anndata

    from sopa.annotation.tangram.annotate import tangram_annotate
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)
    adata_sc = anndata.read(sc_reference_path)

    tangram_annotate(sdata, adata_sc, cell_type_key, bag_size=bag_size)
    sdata.table.write_zarr(Path(sdata_path) / "table" / "table")
