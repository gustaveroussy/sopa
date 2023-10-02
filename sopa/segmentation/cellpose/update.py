import geopandas as gpd
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from ..._constants import SopaKeys
from ..._sdata import get_spatial_image


def add_shapes(sdata: SpatialData, polygons: list[Polygon], image_key: str):
    _, image = get_spatial_image(sdata, image_key)

    geo_df = gpd.GeoDataFrame({"geometry": polygons})
    geo_df.index = image_key + geo_df.index.astype(str)

    geo_df = ShapesModel.parse(geo_df, transformations=get_transformation(image, get_all=True))
    sdata.add_shapes(SopaKeys.CELLPOSE_BOUNDARIES, geo_df, overwrite=True)
