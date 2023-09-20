import typer

app_read = typer.Typer()
option = typer.Option()


@app_read.command()
def qptiff(sdata_path: str, qptiff_path: str = option):
    from sopa.io.qptiff import read_qptiff

    sdata = read_qptiff(qptiff_path)
    sdata.write(sdata_path)


@app_read.command()
def merscope():
    ...
