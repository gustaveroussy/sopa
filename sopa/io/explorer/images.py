import logging
from math import ceil

import numpy as np
import tifffile as tf
import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage, to_multiscale
from spatial_image import SpatialImage

from ._constants import image_metadata

log = logging.getLogger(__name__)


def _astype_uint8(arr: np.ndarray) -> np.ndarray:
    assert np.issubdtype(
        arr.dtype, np.integer
    ), f"The image dtype has to be an integer dtype. Found {arr.dtype}"

    if arr.dtype == np.uint8:
        return arr

    factor = np.iinfo(np.uint8).max / np.iinfo(arr.dtype).max
    return (arr * factor).astype(np.uint8)


def _get_tiles(xarr: xr.DataArray, tile_width: int):
    log.info(f"   Image of shape {xarr.shape}")
    for c in range(xarr.shape[0]):
        for index_x in range(ceil(xarr.shape[2] / tile_width)):
            for index_y in range(ceil(xarr.shape[1] / tile_width)):
                tile = xarr[
                    c,
                    tile_width * index_x : tile_width * (index_x + 1),
                    tile_width * index_y : tile_width * (index_y + 1),
                ].values
                yield _astype_uint8(tile)


def _write_image_level(
    tif: tf.TiffWriter, xarr: xr.DataArray, tile_width, metadata, resolution, **kwargs
):
    tif.write(
        _get_tiles(xarr, tile_width),
        tile=(tile_width, tile_width),
        resolution=(resolution, resolution),
        metadata=metadata,
        shape=xarr.shape,
        dtype=xarr.dtype,
        photometric="minisblack",
        compression="jpeg2000",
        resolutionunit="CENTIMETER",
        **kwargs,
    )


def write_image(
    path: str,
    image: SpatialImage,
    pixelsize: float = 0.2125,
    tile_width: int = 1024,
):
    log.info("Writing multiscale image")

    image_key = image.name
    image: MultiscaleSpatialImage = to_multiscale(image, [2, 2, 2, 2, 2])

    scale_names = list(image.children)
    channel_names = list(map(str, image[scale_names[0]].c.values))

    metadata = image_metadata(channel_names, pixelsize)

    with tf.TiffWriter(path, bigtiff=True) as tif:
        _write_image_level(
            tif,
            image[scale_names[0]][image_key],
            tile_width,
            metadata,
            1e4 / pixelsize,
            subifds=len(scale_names) - 1,
        )

        for i, scale in enumerate(scale_names[1:]):
            _write_image_level(
                tif,
                image[scale][image_key],
                tile_width,
                metadata,
                1e4 * 2 ** (i + 1) / pixelsize,
                subfiletype=1,
            )
