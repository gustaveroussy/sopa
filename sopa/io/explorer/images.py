import logging
from math import ceil

import numpy as np
import tifffile as tf
import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage, to_multiscale
from spatial_image import SpatialImage

from ...utils.image import scale_dtype
from ._constants import image_metadata

log = logging.getLogger(__name__)


class MultiscaleImageWriter:
    photometric = "minisblack"
    compression = "jpeg2000"
    resolutionunit = "CENTIMETER"
    dtype = np.uint8

    def __init__(self, image: MultiscaleSpatialImage, tile_width: int, pixelsize: float):
        self.image = image
        self.tile_width = tile_width
        self.pixelsize = pixelsize

        self.scale_names = list(image.children)
        self.channel_names = list(map(str, image[self.scale_names[0]].c.values))
        self.metadata = image_metadata(self.channel_names, pixelsize)

    def _get_tiles(self, xarr: xr.DataArray):
        for c in range(xarr.shape[0]):
            for index_y in range(ceil(xarr.shape[1] / self.tile_width)):
                for index_x in range(ceil(xarr.shape[2] / self.tile_width)):
                    tile = xarr[
                        c,
                        self.tile_width * index_y : self.tile_width * (index_y + 1),
                        self.tile_width * index_x : self.tile_width * (index_x + 1),
                    ].values
                    yield scale_dtype(tile, self.dtype)

    def _write_image_level(self, tif: tf.TiffWriter, scale_index: int, **kwargs):
        xarr: xr.DataArray = next(iter(self.image[self.scale_names[scale_index]].values()))
        resolution = 1e4 * 2**scale_index / self.pixelsize

        data = self._get_tiles(xarr) if self.lazy else scale_dtype(xarr.values, self.dtype)

        log.info(f"   > Image of shape {xarr.shape}")
        tif.write(
            data,
            tile=(self.tile_width, self.tile_width),
            resolution=(resolution, resolution),
            metadata=self.metadata,
            shape=xarr.shape,
            dtype=self.dtype,
            photometric=self.photometric,
            compression=self.compression,
            resolutionunit=self.resolutionunit,
            **kwargs,
        )

    def __len__(self):
        return len(self.scale_names)

    def write(self, path, lazy=True):
        self.lazy = lazy
        log.info(f"Writing multiscale image (lazy-mode={self.lazy})")

        with tf.TiffWriter(path, bigtiff=True) as tif:
            self._write_image_level(tif, 0, subifds=len(self) - 1)

            for i in range(1, len(self)):
                self._write_image_level(tif, i, subfiletype=1)


def write_image(
    path: str,
    image: SpatialImage | np.ndarray,
    lazy: bool = True,
    tile_width: int = 1024,
    n_subscales: int = 5,
    pixelsize: float = 0.2125,
):
    if isinstance(image, np.ndarray):
        assert len(image.shape) == 3, "Can only write channels with shape (C,Y,X)"
        log.info(f"Converting image of shape {image.shape} into a SpatialImage (with dims: C,Y,X)")
        image = SpatialImage(image, dims=["c", "y", "x"], name="image")

    image: MultiscaleSpatialImage = to_multiscale(image, [2] * n_subscales)

    image_writer = MultiscaleImageWriter(image, pixelsize=pixelsize, tile_width=tile_width)
    image_writer.write(path, lazy=lazy)
