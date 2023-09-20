import geopandas as gpd
import pandas as pd
from anndata import AnnData
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel, TableModel

from ..utils.utils import _get_spatial_image
from .shapes import average


def update(sdata: SpatialData, polygons: list[Polygon], image_key: str):
    image = _get_spatial_image(sdata, image_key)

    geo_df = gpd.GeoDataFrame(
        {
            "geometry": polygons,
            "x": [poly.centroid.x for poly in polygons],
            "y": [poly.centroid.y for poly in polygons],
        }
    )
    geo_df.index = image_key + geo_df.index.astype(str)

    geo_df = ShapesModel.parse(geo_df, transformations=image.transform)
    sdata.add_shapes("polygons", geo_df)

    mean_intensities = average(image, polygons)
    adata = AnnData(
        mean_intensities,
        dtype=mean_intensities.dtype,
        var=pd.DataFrame(index=image.c),
        obs=pd.DataFrame(index=geo_df.index),
    )

    adata.obsm["spatial"] = geo_df[["x", "y"]].values
    adata.obs["region"] = pd.Series("polygons", index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(image_key, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(image_key, index=adata.obs_names, dtype="category")
    adata.obs["cell_id"] = geo_df.index

    adata = TableModel.parse(
        adata,
        region_key="region",
        region="polygons",
        instance_key="cell_id",
    )

    sdata.table = adata
