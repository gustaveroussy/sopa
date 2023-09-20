import argparse
import re
import shutil
from pathlib import Path

import dask.array as da
import tifffile as tf
import toml
import xarray as xr
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity


def get_channel_name(description):
    return re.search(r"<Name>(.*?)</Name>", description).group(1)


def read_qptiff(
    path: Path, channels_renaming: dict | None = None, image_models_kwargs: dict | None = None
) -> SpatialData:
    image_models_kwargs = {} if image_models_kwargs is None else image_models_kwargs
    if "chunks" not in image_models_kwargs:
        image_models_kwargs["chunks"] = (1, 4096, 4096)

    with tf.TiffFile(path) as tif:
        page_series = tif.series[0]
        names = [get_channel_name(page.description) for page in page_series]

        if channels_renaming is not None:
            names = [channels_renaming.get(name, name) for name in names]

        image_name = Path(path).absolute().stem
        image = Image2DModel.parse(
            da.from_array(page_series.asarray(), chunks=image_models_kwargs["chunks"]),
            dims=list(page_series._axes.lower()),
            transformations={"pixels": Identity()},
            c_coords=names,
            **image_models_kwargs,
        )

        return SpatialData(images={image_name: image})


def main(args):
    path = Path(args.path)

    config = toml.load(args.config)

    sdata = read_qptiff(path, channels_renaming=config["reader"]["channels_renaming"])

    print(sdata)
    sdata.write(args.sdata_path)


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
        "-s",
        "--sdata_path",
        type=str,
        required=True,
        help="Path to the zarr sdata",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="Path to the config file",
    )

    main(parser.parse_args())
