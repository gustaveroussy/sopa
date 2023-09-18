import click
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata
from matplotlib.widgets import PolygonSelector
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel

from .._constants import ROI
from .image import resize
from .utils import _get_spatial_image

HELPER = """Enclose cells within a polygon. Helper:
    - Click on the plot to add a polygon vertex
    - Press the 'esc' key to start a new polygon
    - Try holding the 'ctrl' key to move a single vertex
    - Once the polygon is finished and overlaid in red, you can close the window
"""

VALID_N_CHANNELS = [1, 3]


class _Selector:
    def __init__(self, ax):
        self.poly = PolygonSelector(ax, self.onselect, draw_bounding_box=True)
        print(HELPER)
        plt.show()

    def onselect(self, vertices):
        self.vertices = np.array(vertices)
        print(f"Selected polygon with {len(self.vertices)} vertices.")

    def disconnect(self):
        self.poly.disconnect_events()


def xarr_selector(
    sdata: SpatialData,
    channels: list[str],
    scale_factor: float = 10,
    margin_ratio: float = 0.1,
):
    _, spatial_image = _get_spatial_image(sdata)
    image = spatial_image.transpose("y", "x", "c")

    if len(channels):
        assert (
            len(channels) in VALID_N_CHANNELS
        ), f"Number of channels provided must be in: {', '.join(VALID_N_CHANNELS)}"
        image = image.sel(c=channels)
    else:
        assert (
            len(image.coords["c"]) in VALID_N_CHANNELS
        ), f"Choose one or three channels among {image.c.values} by using the -c argument"

    image = resize(image, scale_factor).compute()

    _, ax = plt.subplots()
    ax.imshow(image)

    dy, dx, *_ = image.shape
    plt.xlim(-margin_ratio * dx, dx + margin_ratio * dx)
    plt.ylim(dy + margin_ratio * dy, -margin_ratio * dy)

    selector = _Selector(ax)

    polygon = Polygon(selector.vertices * scale_factor)

    geo_df = gpd.GeoDataFrame({"geometry": [polygon]})
    geo_df = ShapesModel.parse(geo_df, transformations=spatial_image.transform)
    sdata.add_shapes(ROI, geo_df)


# def cells_selector(metadata_path: str, output_path: str, x: str = "center_x", y: str = "center_y"):
#     df = pd.read_csv(metadata_path)

#     _, ax = plt.subplots()
#     ax.scatter(
#         df[x],
#         df[y],
#         marker=".",
#         rasterized=True,
#         s=0.05,
#     )

#     selector = _Selector(ax)

#     np.savetxt(output_path, selector.vertices)


@click.command()
@click.option(
    "-i",
    "--input_path",
    type=str,
    help="Path to input sdata.zarr directory",
)
@click.option(
    "-c",
    "--channels",
    type=str,
    multiple=True,
    default=None,
    help="List of channels name to be displayed",
)
def main(input_path, channels):
    sdata = spatialdata.read_zarr(input_path)

    xarr_selector(sdata, list(channels))


if __name__ == "__main__":
    main()
