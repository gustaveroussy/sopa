import argparse
import re
import shutil
from pathlib import Path
from typing import List

import tifffile as tf

from .xenium import write_ome_tif


def get_channel_name(description):
    return re.search(r"<Name>(.*?)</Name>", description).group(1)


def read_series(path: Path) -> List[tf.TiffPageSeries]:
    with tf.TiffFile(path) as tif:
        return list(reversed(sorted(tif.series[0], key=lambda p: p.size)))


def write_zarr(
    path: Path,
    series: List[tf.TiffPageSeries],
    names: List[str],
    overwrite: bool = True,
) -> None:
    import dask.array as da
    import xarray as xr

    dask_array = da.asarray(series[0].asarray())
    xarr = xr.DataArray(
        dask_array, dims=list(series[0]._axes.lower()), coords={"c": names}
    )
    ds = xr.Dataset({"image": xarr})

    if path.exists():
        assert overwrite, f"Path {path} exists and overwrite is False"
        shutil.rmtree(path)

    print("Saving xarray")
    ds.to_zarr(path)


def main(args):
    path, output = Path(args.path), Path(args.output)

    assert not output.exists(), f"Output path {output} already exists"

    series = read_series(path)
    names = [get_channel_name(page.description) for page in series[0]._pages]

    write_ome_tif(output, series, names)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        required=True,
        help="Path to the qptiff file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the morphology.ome.tif file",
    )

    main(parser.parse_args())
