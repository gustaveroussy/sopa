import cv2
import numpy as np
import geopandas as gpd
from sopa.io import read_wsi
from xarray import DataArray
from spatialdata import SpatialData
from shapely.geometry import Polygon
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

def hsv_otsu(
    sdata: SpatialData,
    lvl: int=-1,
    blur_k: int=5, 
    open_k: int=5, 
    close_k: int=5, 
    drop_threshold: int=0.01
):

    level_keys = list(sdata.images.keys())
    image = sdata[level_keys[lvl]]

    thumbnail = np.array(image.transpose('y','x','c'))
    thumbnail_hsv = cv2.cvtColor(thumbnail, cv2.COLOR_RGB2HSV)
    thumbnail_hsv_blurred = cv2.medianBlur(thumbnail_hsv[:, :, 1], blur_k)
    _, mask = cv2.threshold(
        thumbnail_hsv_blurred, 0, 255, cv2.THRESH_OTSU + cv2.THRESH_BINARY
    )

    mask_open = cv2.morphologyEx(
        mask, cv2.MORPH_OPEN, np.ones((open_k, open_k), np.uint8)
    )
    mask_open_close = cv2.morphologyEx(
        mask_open, cv2.MORPH_CLOSE, np.ones((close_k, close_k), np.uint8)
    )

    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(
        mask_open_close, 4, cv2.CV_32S
    )

    contours = []
    for i in range(1, num_labels):
        if stats[i, 4] > drop_threshold * np.prod(mask_open_close.shape):
            cc = cv2.findContours(
                np.array(labels == i, dtype="uint8"),
                cv2.RETR_TREE,
                cv2.CHAIN_APPROX_NONE,
            )[0][0]
            c_closed = np.array(list(cc) + [cc[0]])
            contours.extend([c_closed.squeeze()])

    polygons = [Polygon(contour) for contour in contours]
    geo_df = gpd.GeoDataFrame(geometry=polygons)

    geo_df = ShapesModel.parse(
        geo_df,
        transformations=get_transformation(image, get_all=True)
    )
    sdata.add_shapes('tissue', geo_df)

if __name__ == '__main__':
    img = read_wsi('CMU-1.ndpi')
    hsv_otsu(img)
    
    # import spatialdata_plot
    # img\
    #     .pl.render_images('4')\
    #     .pl.render_shapes(outline_color='green',fill_alpha=0, outline=True)\
    #     .pl.show(save='CMU-1-tissue_segmentation.png')
