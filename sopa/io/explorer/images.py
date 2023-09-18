import numpy as np
import tifffile as tf
from multiscale_spatial_image import MultiscaleSpatialImage

from ._constants import image_metadata, image_options


def _astype_uint8(arr: np.ndarray) -> np.ndarray:
    print(f"   Image of shape {arr.shape}")
    assert np.issubdtype(
        arr.dtype, np.integer
    ), f"The image dtype has to be an integer dtype. Found {arr.dtype}"

    if arr.dtype == np.uint8:
        return arr

    factor = np.iinfo(np.uint8).max / np.iinfo(arr.dtype).max
    return (arr * factor).astype(np.uint8)


def write_multiscale(
    path: str,
    multiscale: MultiscaleSpatialImage,
    pixelsize: float = 0.2125,
):
    print("Writing multiscale image")

    assert isinstance(
        multiscale, MultiscaleSpatialImage
    ), f"For now, only `MultiscaleSpatialImage` is supported. Found {type(multiscale)}."

    scale_names = list(multiscale.children)
    channel_names = list(multiscale[scale_names[0]].c.values)

    metadata = image_metadata(channel_names, pixelsize)

    # TODO : make it memory efficient
    with tf.TiffWriter(path, bigtiff=True) as tif:
        tif.write(
            _astype_uint8(multiscale[scale_names[0]]["image"].values),
            subifds=len(scale_names) - 1,
            resolution=(1e4 / pixelsize, 1e4 / pixelsize),
            metadata=metadata,
            **image_options(),
        )

        for i, scale in enumerate(scale_names[1:]):
            tif.write(
                _astype_uint8(multiscale[scale]["image"].values),
                subfiletype=1,
                resolution=(
                    1e4 * 2 ** (i + 1) / pixelsize,
                    1e4 * 2 ** (i + 1) / pixelsize,
                ),
                **image_options(),
            )
