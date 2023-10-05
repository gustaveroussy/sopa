import logging
from math import ceil

import numpy as np
import tifffile as tf
import xarray as xr
from multiscale_spatial_image import to_multiscale
from spatial_image import SpatialImage

from ._constants import image_metadata, image_options

log = logging.getLogger(__name__)


def _astype_uint8(arr: np.ndarray) -> np.ndarray:
    log.info(f"   Image of shape {arr.shape}")
    assert np.issubdtype(
        arr.dtype, np.integer
    ), f"The image dtype has to be an integer dtype. Found {arr.dtype}"

    if arr.dtype == np.uint8:
        return arr

    factor = np.iinfo(np.uint8).max / np.iinfo(arr.dtype).max
    return (arr * factor).astype(np.uint8)


def write_image(
    path: str,
    image: SpatialImage,
    image_key: str,
    pixelsize: float = 0.2125,
):
    log.info("Writing multiscale image")

    image = to_multiscale(image, [2, 2, 2, 2, 2])

    scale_names = list(image.children)
    channel_names = list(map(str, image[scale_names[0]].c.values))

    metadata = image_metadata(channel_names, pixelsize)

    # TODO : switch to lazy
    with tf.TiffWriter(path, bigtiff=True) as tif:
        tif.write(
            _astype_uint8(image[scale_names[0]][image_key].values),
            subifds=len(scale_names) - 1,
            resolution=(1e4 / pixelsize, 1e4 / pixelsize),
            metadata=metadata,
            **image_options(),
        )

        for i, scale in enumerate(scale_names[1:]):
            tif.write(
                _astype_uint8(image[scale][image_key].values),
                subfiletype=1,
                resolution=(
                    1e4 * 2 ** (i + 1) / pixelsize,
                    1e4 * 2 ** (i + 1) / pixelsize,
                ),
                **image_options(),
            )


def get_tiles(image: xr.DataArray, tile_width: int):
    # TODO: WIP: do it with dask
    for c in range(image.shape[0]):
        for x in range(ceil(image.shape[2] / tile_width)):
            for y in range(ceil(image.shape[1] / tile_width)):
                tile = image[
                    c, tile_width * x : tile_width * (x + 1), tile_width * y : tile_width * (y + 1)
                ]
                yield _astype_uint8(tile)


def _write_image_level(tif: tf.TiffWriter, image: xr.DataArray, resolution, metadata, **kwargs):
    tif.write(
        get_tiles(image),
        resolution=(resolution, resolution),
        metadata=metadata,
        shape=image.shape,
        dtype=image.dtype,
        **image_options(),
        **kwargs,
    )


def write_image_lazy(
    path: str,
    image: SpatialImage,
    image_key: str,
    pixelsize: float = 0.2125,
):
    log.info("Writing multiscale image")

    image = to_multiscale(image, [2, 2, 2, 2, 2])

    scale_names = list(image.children)
    channel_names = list(map(str, image[scale_names[0]].c.values))

    metadata = image_metadata(channel_names, pixelsize)

    with tf.TiffWriter(path, bigtiff=True) as tif:
        _write_image_level(
            tif,
            image[scale_names[0]][image_key],
            1e4 / pixelsize,
            metadata,
            subifds=len(scale_names) - 1,
        )

        for i, scale in enumerate(scale_names[1:]):
            _write_image_level(
                tif,
                image[scale][image_key],
                1e4 * 2 ** (i + 1) / pixelsize,
                metadata,
                subfiletype=1,
            )
