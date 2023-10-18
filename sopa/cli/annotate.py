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

    import spatialdata

    from sopa.annotation.fluorescence import higher_z_score

    sdata = spatialdata.read_zarr(sdata_path)

    assert sdata.table is not None, f"Annotation requires `sdata.table` to be not None"

    higher_z_score(sdata.table, marker_cell_dict, key)
    sdata.table.write_zarr(Path(sdata_path) / "table" / "table")


@app_annotate.command()
def tangram(
    sdata_path: str,
    reference_path: str,
    key: str = "cell_type",
):
    ...  # TODO
