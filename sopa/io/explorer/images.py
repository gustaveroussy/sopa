import logging
import re
from math import ceil

import numpy as np
import tifffile as tf
import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage, to_multiscale
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Affine
from tqdm import tqdm

from ..._sdata import get_intrinsic_cs, get_spatial_image
from ...utils.image import resize_numpy, scale_dtype
from ..imaging import _default_image_models_kwargs
from ._constants import ExplorerConstants, FileNames, image_metadata
from .utils import explorer_file_path

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
        self.channel_names = _set_colors(self.channel_names)
        self.metadata = image_metadata(self.channel_names, pixelsize)
        self.data = None

        self.lazy = True
        self.ram_threshold_gb = None

    def _n_tiles_axis(self, xarr: xr.DataArray, axis: int) -> int:
        return ceil(xarr.shape[axis] / self.tile_width)

    def _get_tiles(self, xarr: xr.DataArray):
        for c in range(xarr.shape[0]):
            for index_y in range(self._n_tiles_axis(xarr, 1)):
                for index_x in range(self._n_tiles_axis(xarr, 2)):
                    tile = xarr[
                        c,
                        self.tile_width * index_y : self.tile_width * (index_y + 1),
                        self.tile_width * index_x : self.tile_width * (index_x + 1),
                    ].values
                    yield self._scale(tile)

    def _should_load_memory(self, shape: tuple[int, int, int], dtype: np.dtype):
        if not self.lazy:
            return True

        if self.ram_threshold_gb is None:
            return False

        itemsize = max(np.dtype(dtype).itemsize, np.dtype(self.dtype).itemsize)
        size = shape[0] * shape[1] * shape[2] * itemsize

        return size <= self.ram_threshold_gb * 1024**3

    def _scale(self, array: np.ndarray):
        return scale_dtype(array, self.dtype).clip(0, np.iinfo(np.int8).max)

    def _write_image_level(self, tif: tf.TiffWriter, scale_index: int, **kwargs):
        xarr: xr.DataArray = next(iter(self.image[self.scale_names[scale_index]].values()))
        resolution = 1e4 * 2**scale_index / self.pixelsize

        if not self._should_load_memory(xarr.shape, xarr.dtype):
            n_tiles = xarr.shape[0] * self._n_tiles_axis(xarr, 1) * self._n_tiles_axis(xarr, 2)
            data = self._get_tiles(xarr)
            data = iter(tqdm(data, total=n_tiles - 1, desc="Writing tiles"))
        else:
            if self.data is not None:
                self.data = resize_numpy(self.data, 2, xarr.dims, xarr.shape)
            else:
                log.info(f"   (Loading image of shape {xarr.shape}) in memory")
                self.data = self._scale(xarr.values)

            data = self.data

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

    def procedure(self):
        if not self.lazy:
            return "in-memory (consider lazy procedure if it crashes because of RAM)"
        if self.ram_threshold_gb is None:
            return "lazy (slower but low RAM usage)"
        return "semi-lazy (load in memory when possible)"

    def write(self, path, lazy=True, ram_threshold_gb=None):
        self.lazy = lazy
        self.ram_threshold_gb = ram_threshold_gb

        log.info(f"Writing multiscale image with procedure={self.procedure()}")

        with tf.TiffWriter(path, bigtiff=True) as tif:
            self._write_image_level(tif, 0, subifds=len(self) - 1)

            for i in range(1, len(self)):
                self._write_image_level(tif, i, subfiletype=1)


def _to_color(channel_name: str, is_wavelength: bool, colors_iterator: list):
    if is_wavelength:
        return channel_name
    if channel_name in ExplorerConstants.KNOWN_CHANNELS:
        return f"{channel_name} (color={ExplorerConstants.KNOWN_CHANNELS[channel_name]})"
    return f"{channel_name} (color={colors_iterator.pop()})"


def _set_colors(channel_names: list[str]) -> list[str]:
    """
    Trick to provide a color to all channels on the Xenium Explorer.

    But some channels colors are set to white by default. This functions allows to color these
    channels with an available wavelength color (e.g., `550`).
    """
    existing_wavelength = [
        bool(re.search(r"(?<![0-9])[0-9]{3}(?![0-9])", c)) for c in channel_names
    ]
    valid_colors = [c for c in ExplorerConstants.COLORS if c != ExplorerConstants.NUCLEUS_COLOR]
    n_missing = sum(
        not is_wavelength and not c in ExplorerConstants.KNOWN_CHANNELS
        for c, is_wavelength in zip(channel_names, existing_wavelength)
    )
    colors_iterator: list = np.repeat(valid_colors, ceil(n_missing / len(valid_colors))).tolist()

    return [
        _to_color(c, is_wavelength, colors_iterator)
        for c, is_wavelength in zip(channel_names, existing_wavelength)
    ]


def write_image(
    path: str,
    image: SpatialImage | np.ndarray,
    lazy: bool = True,
    tile_width: int = 1024,
    n_subscales: int = 5,
    pixelsize: float = 0.2125,
    ram_threshold_gb: int | None = 4,
    is_dir: bool = True,
):
    """Convert an image into a `morphology.ome.tif` file that can be read by the Xenium Explorer

    Args:
        path: Path to the Xenium Explorer directory where the image will be written
        image: Image of shape `(C, Y, X)`
        lazy: If `False`, the image will not be read in-memory (except if the image size is below `ram_threshold_gb`). If `True`, all the images levels are always loaded in-memory.
        tile_width: Xenium tile width (do not update).
        n_subscales: Number of sub-scales in the pyramidal image.
        pixelsize: Xenium pixel size (do not update).
        ram_threshold_gb: If an image (of any level of the pyramid) is below this threshold, it will be loaded in-memory.
        is_dir: If `False`, then `path` is a path to a single file, not to the Xenium Explorer directory.
    """
    path = explorer_file_path(path, FileNames.IMAGE, is_dir)

    if isinstance(image, np.ndarray):
        assert len(image.shape) == 3, "Can only write channels with shape (C,Y,X)"
        log.info(f"Converting image of shape {image.shape} into a SpatialImage (with dims: C,Y,X)")
        image = SpatialImage(image, dims=["c", "y", "x"], name="image")

    image: MultiscaleSpatialImage = to_multiscale(image, [2] * n_subscales)

    image_writer = MultiscaleImageWriter(image, pixelsize=pixelsize, tile_width=tile_width)
    image_writer.write(path, lazy=lazy, ram_threshold_gb=ram_threshold_gb)


def align(
    sdata: SpatialData,
    image: SpatialImage,
    transformation_matrix_path: str,
    image_key: str = None,
    image_models_kwargs: dict | None = None,
    overwrite: bool = False,
):
    """Add an image to the `SpatialData` object after alignment with the Xenium Explorer.

    Args:
        sdata: A `SpatialData` object
        image: A `SpatialImage` object. Note that `image.name` is used as the key for the aligned image.
        transformation_matrix_path: Path to the `.csv` transformation matrix exported from the Xenium Explorer
        image_key: Optional name of the image on which it has been aligned. Required if multiple images in the `SpatialData` object.
        image_models_kwargs: Kwargs to the `Image2DModel` model.
        overwrite: Whether to overwrite the image, if already existing.
    """
    image_models_kwargs = _default_image_models_kwargs(image_models_kwargs)

    to_pixel = Affine(
        np.genfromtxt(transformation_matrix_path, delimiter=","),
        input_axes=("x", "y"),
        output_axes=("x", "y"),
    )

    default_image = get_spatial_image(sdata, image_key)
    pixel_cs = get_intrinsic_cs(sdata, default_image)

    image = Image2DModel.parse(
        image,
        dims=("c", "y", "x"),
        transformations={pixel_cs: to_pixel},
        c_coords=image.coords["c"].values,
        **image_models_kwargs,
    )

    log.info(f"Adding image {image.name}:\n{image}")
    sdata.add_image(image.name, image, overwrite=overwrite)
