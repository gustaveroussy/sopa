import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel, TableModel

from ..utils.utils import _get_spatial_image
from .shapes import average


def update(sdata: SpatialData, polygons: list[Polygon], image_key: str):
    _, image = _get_spatial_image(sdata, image_key)

    geo_df = gpd.GeoDataFrame({"geometry": polygons})
    geo_df.index = image_key + geo_df.index.astype(str)

    geo_df = ShapesModel.parse(geo_df, transformations=image.transform)
    sdata.add_shapes("polygons", geo_df)  # TODO: not hardcoded name


def aggregate(sdata: SpatialData, table: AnnData | None, intensity_mean: bool):
    geo_df = sdata["polygons"]

    image_key, image = _get_spatial_image(sdata)

    assert (
        intensity_mean or table is not None
    ), f"You must choose at least one aggregation: transcripts or fluorescence intensities"

    if intensity_mean:
        mean_intensities = average(image, geo_df.geometry)

    if table is None:
        table = AnnData(
            mean_intensities,
            dtype=mean_intensities.dtype,
            var=pd.DataFrame(index=image.c),
            obs=pd.DataFrame(index=geo_df.index),
        )
    else:
        table.obsm["intensities"] = pd.DataFrame(mean_intensities, columns=list(image.c))

    centroids = geo_df.centroid

    table.obsm["spatial"] = np.array([[centroid.x, centroid.y] for centroid in centroids])
    table.obs["region"] = pd.Series("polygons", index=table.obs_names, dtype="category")
    table.obs["slide"] = pd.Series(image_key, index=table.obs_names, dtype="category")
    table.obs["dataset_id"] = pd.Series(image_key, index=table.obs_names, dtype="category")
    table.obs["cell_id"] = geo_df.index

    table = TableModel.parse(
        table,
        region_key="region",
        region="polygons",
        instance_key="cell_id",
    )

    sdata.table = table
