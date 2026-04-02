import geopandas as gpd
import pandas as pd
from anndata import AnnData
from spatialdata import SpatialData
from spatialdata.models import TableModel

from ..constants import ATTRS_KEY, SopaKeys
from ..io.explorer.utils import str_cell_id
from ..utils import add_spatial_element


def add_parsed_table(
    sdata: SpatialData,
    table: AnnData,
    geo_df: gpd.GeoDataFrame,
    shapes_key: str,
    table_key: str,
    image_key: str | None = None,
    add_shapes: bool = True,
):
    table = parse_table(table, geo_df, shapes_key, image_key)

    if add_shapes:
        add_spatial_element(sdata, shapes_key, geo_df)

    add_spatial_element(sdata, table_key, table)


def parse_table(table: AnnData, geo_df: gpd.GeoDataFrame, shapes_key: str, image_key: str | None = None) -> AnnData:
    table.obs_names = list(map(str_cell_id, range(table.n_obs)))
    geo_df.index = table.obs_names.to_list()

    table.obsm["spatial"] = geo_df.centroid.get_coordinates().values
    table.obs[SopaKeys.REGION_KEY] = pd.Series(shapes_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.SLIDE_KEY] = pd.Series(image_key or "None", index=table.obs_names, dtype="category")
    table.obs[SopaKeys.INSTANCE_KEY] = geo_df.index
    table.obs[SopaKeys.AREA_OBS] = geo_df.area.values

    if ATTRS_KEY in table.uns:
        del table.uns[ATTRS_KEY]

    table = TableModel.parse(
        table,
        region_key=SopaKeys.REGION_KEY,
        region=shapes_key,
        instance_key=SopaKeys.INSTANCE_KEY,
    )

    return table
