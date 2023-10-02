import ast

import typer

app_read = typer.Typer()
option = typer.Option()


@app_read.command()
def qptiff(
    sdata_path: str,
    qptiff_path: str,
    channels_renaming: str = typer.Option(default={}, callback=ast.literal_eval),
):
    from sopa.io.qptiff import read_qptiff

    sdata = read_qptiff(qptiff_path, channels_renaming=channels_renaming)
    sdata.write(sdata_path)


@app_read.command()
def merscope(path: str, sdata_path: str = None):
    from pathlib import Path

    from spatialdata import SpatialData

    from sopa import io

    path = Path(path)

    sdata_path = path.with_suffix(".zarr") if sdata_path is None else sdata_path

    sdata: SpatialData = io.merscope(path)
    sdata.write(sdata_path)
