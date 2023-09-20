import ast

import typer

app_read = typer.Typer()
option = typer.Option()


@app_read.command()
def qptiff(
    sdata_path: str,
    qptiff_path: str,
    channels_renaming: str = typer.Option(default=None, callback=ast.literal_eval),
):
    from sopa.io.qptiff import read_qptiff

    sdata = read_qptiff(qptiff_path, channels_renaming=channels_renaming)
    sdata.write(sdata_path)


@app_read.command()
def merscope():
    ...
