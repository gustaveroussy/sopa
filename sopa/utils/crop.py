### Deprecated functionnalities

import logging

import dask.array as da
import dask_image.ndinterp
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import zarr
from matplotlib.widgets import PolygonSelector
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from xarray import DataArray

from .._constants import SopaKeys
from .utils import get_spatial_image

log = logging.getLogger(__name__)

HELPER = """Enclose cells within a polygon. Helper:
    - Click on the plot to add a polygon vertex
    - Press the 'esc' key to start a new polygon
    - Try holding the 'ctrl' key to move a single vertex
    - Once the polygon is finished and overlaid in red, you can close the window
"""

VALID_N_CHANNELS = [1, 3]


def _prepare(sdata: SpatialData, channels: list[str], scale_factor: float):
    image_key, spatial_image = get_spatial_image(sdata, return_key=True)
    image = spatial_image.transpose("y", "x", "c")

    if channels is not None and len(channels):
        assert (
            len(channels) in VALID_N_CHANNELS
        ), f"Number of channels provided must be in: {', '.join(VALID_N_CHANNELS)}"
        image = image.sel(c=channels)
    else:
        assert (
            len(image.coords["c"]) in VALID_N_CHANNELS
        ), f"Choose one or three channels among {image.c.values} by using the --channels argument"

    log.info(f"Resizing image by a factor of {scale_factor}")
    return image_key, _resize_dataarray(image, scale_factor).compute()


class _Selector:
    def __init__(self, ax):
        self.poly = PolygonSelector(ax, self.onselect, draw_bounding_box=True)
        log.info(HELPER)
        plt.show()

    def onselect(self, vertices):
        self.vertices = np.array(vertices)
        log.info(f"Selected polygon with {len(self.vertices)} vertices.")

    def disconnect(self):
        self.poly.disconnect_events()


def _draw_polygon(image: np.ndarray, scale_factor: float, margin_ratio: float):
    _, ax = plt.subplots()
    ax.imshow(image)

    dy, dx, *_ = image.shape
    plt.xlim(-margin_ratio * dx, dx + margin_ratio * dx)
    plt.ylim(dy + margin_ratio * dy, -margin_ratio * dy)

    selector = _Selector(ax)

    return Polygon(selector.vertices * scale_factor)


def intermediate_selection(intermediate_image: str, intermediate_polygon: str, margin_ratio: float = 0.1):
    log.info(f"Reading intermediate image {intermediate_image}")

    z = zarr.open(intermediate_image, mode="r")
    image = z["image"][:]

    polygon = _draw_polygon(image, z.attrs["scale_factor"], margin_ratio)

    with zarr.ZipStore(intermediate_polygon, mode="w") as store:
        g = zarr.group(store=store)
        g.attrs.put({"image_key": z.attrs["image_key"], "element_type": "images"})

        coords = np.array(polygon.exterior.coords)
        g.array("polygon", coords, dtype=coords.dtype, chunks=coords.shape)


def polygon_selection(
    sdata: SpatialData,
    intermediate_image: str | None = None,
    intermediate_polygon: str | None = None,
    channels: list[str] | None = None,
    scale_factor: float = 10,
    margin_ratio: float = 0.1,
):
    """Crop an image based on a user-defined polygon (interactive mode).

    Warning:
        This function is deprecated. Using `napari-spatialdata` instead.
        Provide the name `"region_of_interest"` to your selected ROI.

    Args:
        sdata: A `SpatialData` object
        intermediate_image: Path to the intermediate image, with a `.zip` extension. Use this only if the interactive mode is not available
        intermediate_polygon: Path to the intermediate polygon, with a `.zip` extension. Use this locally, after downloading the `intermediate_image`
        channels: List of channel names to be displayed. Optional if there are already only 1 or 3 channels.
        scale_factor: Resize the image by this value (high value for a lower memory usage)
        margin_ratio: Ratio of the image margin on the display (compared to the image size)
    """
    log.warning("This function is deprecated. Use `napari-spatialdata` instead.")

    if intermediate_polygon is None:
        image_key, image = _prepare(sdata, channels, scale_factor)

        if intermediate_image is not None:
            log.info(f"Resized image will be saved to {intermediate_image}")
            with zarr.ZipStore(intermediate_image, mode="w") as store:
                g = zarr.group(store=store)
                g.attrs.put({"scale_factor": scale_factor, "image_key": image_key})
                g.array("image", image, dtype=image.dtype, chunks=image.shape)
            return

        polygon = _draw_polygon(image, scale_factor, margin_ratio)
    else:
        log.info(f"Reading polygon at path {intermediate_polygon}")
        z = zarr.open(intermediate_polygon, mode="r")
        polygon = Polygon(z["polygon"][:])
        image_key = z.attrs["image_key"]

        image = get_spatial_image(sdata, image_key)

    geo_df = gpd.GeoDataFrame(geometry=[polygon])

    geo_df = ShapesModel.parse(geo_df, transformations=get_transformation(sdata[image_key], get_all=True).copy())
    sdata.shapes[SopaKeys.ROI] = geo_df
    if sdata.is_backed():
        sdata.write_element(SopaKeys.ROI, overwrite=True)

    log.info(f"Polygon saved in sdata['{SopaKeys.ROI}']")


def _resize_dataarray(xarr: DataArray, scale_factor: float) -> da.Array:
    """Resize a xarray image

    Args:
        xarr: A `xarray` array
        scale_factor: Scale factor of resizing, e.g. `2` will decrease the width by 2

    Returns:
        Resized dask array
    """
    resize_dims = [dim in ["x", "y"] for dim in xarr.dims]
    transform = np.diag([scale_factor if resize_dim else 1 for resize_dim in resize_dims])
    output_shape = [size // scale_factor if resize_dim else size for size, resize_dim in zip(xarr.shape, resize_dims)]

    return dask_image.ndinterp.affine_transform(xarr.data, matrix=transform, output_shape=output_shape)
