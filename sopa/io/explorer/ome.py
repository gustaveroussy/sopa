from pathlib import Path

import cv2
import numpy as np
import tifffile as tf
from multiscale_spatial_image import MultiscaleSpatialImage

OPTIONS = dict(
    photometric="minisblack",
    tile=(1024, 1024),
    compression="jpeg2000",
    resolutionunit="CENTIMETER",
)


def get_metadata(channel_names, pixelsize) -> dict:
    return {
        "SignificantBits": 8,
        "PhysicalSizeX": pixelsize,
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": pixelsize,
        "PhysicalSizeYUnit": "µm",
        "Channel": {"Name": channel_names},
    }


def img_resize(img, scale_factor):
    width = int(np.floor(img.shape[1] * scale_factor))
    height = int(np.floor(img.shape[0] * scale_factor))
    return cv2.resize(img, (width, height), interpolation=cv2.INTER_AREA)


def write_np_ome_tif(
    path,
    image: np.ndarray,
    subresolutions: int = 7,
    channel_names: list[str] = None,
    pixelsize: float = 0.2125,
):
    if channel_names is None:
        channel_names = list(range(image.shape[0]))

    metadata = get_metadata(channel_names, pixelsize)

    with tf.TiffWriter(path, bigtiff=True) as tif:
        tif.write(
            np.moveaxis(image, -1, 0),
            subifds=subresolutions,
            resolution=(1e4 / pixelsize, 1e4 / pixelsize),
            metadata=metadata,
            **OPTIONS,
        )

        for i in range(1, subresolutions + 1):
            scale = 1 / 2**i
            downsample = img_resize(image, scale)
            tif.write(
                np.moveaxis(downsample, -1, 0),
                subfiletype=1,
                resolution=(1e4 / scale / pixelsize, 1e4 / scale / pixelsize),
                **OPTIONS,
            )


def to_uint8(arr):
    return (arr // 256).astype(np.uint8)


def write_multiscale(
    output_path: Path,
    multiscale: MultiscaleSpatialImage,
    pixelsize: float = 0.2125,
):
    scale_names = list(multiscale.children)
    channel_names = list(multiscale[scale_names[0]].c.values)

    metadata = get_metadata(channel_names, pixelsize)

    with tf.TiffWriter(output_path, bigtiff=True) as tif:
        tif.write(
            to_uint8(multiscale[scale_names[0]]["image"].values),
            subifds=len(scale_names) - 1,
            resolution=(1e4 / pixelsize, 1e4 / pixelsize),
            metadata=metadata,
            **OPTIONS,
        )

        for i, scale in enumerate(scale_names[1:]):
            tif.write(
                to_uint8(multiscale[scale]["image"].values),
                subfiletype=1,
                resolution=(
                    1e4 * 2 ** (i + 1) / pixelsize,
                    1e4 * 2 ** (i + 1) / pixelsize,
                ),
                **OPTIONS,
            )


def write_ome_tif(
    output_path: Path,
    series: list[tf.TiffPageSeries],
    channel_names: list[str] = None,
    pixelsize: float = 0.2125,
) -> None:
    if channel_names is None:
        channel_names = list(range(series[0].shape[0]))

    metadata = get_metadata(channel_names, pixelsize)

    with tf.TiffWriter(output_path, bigtiff=True) as tif:
        tif.write(
            series[0].asarray(),
            subifds=len(series) - 1,
            resolution=(1e4 / pixelsize, 1e4 / pixelsize),
            metadata=metadata,
            **OPTIONS,
        )

        for i in range(1, len(series)):
            print(f"Writing sub-resolution {i}...")
            tif.write(
                series[i].asarray(),
                subfiletype=1,
                resolution=(1e4 * 2**i / pixelsize, 1e4 * 2**i / pixelsize),
                **OPTIONS,
            )
