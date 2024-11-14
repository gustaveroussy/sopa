from functools import partial

import dask.dataframe as dd
import geopandas as gpd
import pandas as pd
from spatialdata import SpatialData

from .._constants import SopaAttrs
from ..utils import get_boundaries, get_spatial_element, to_intrinsic


def sjoin(
    sdata: SpatialData,
    left_element: str | gpd.GeoDataFrame,
    right_element: str | gpd.GeoDataFrame,
    how: str = "left",
    target_coordinate_system: str | None = None,
    **kwargs: int,
) -> gpd.GeoDataFrame:
    """Spatial join of two `shapes` GeoDataFrames, as in [geopandas.sjoin](https://geopandas.org/en/stable/docs/reference/api/geopandas.sjoin.html).

    Shapes are automatically aligned on the same coordinate system (which can be chosen using the `target_coordinate_system` argument).

    Args:
        sdata: A `SpatialData` object
        left_element: The name of a GeoDataFrame in `sdata`, or the GeoDataFrame itself
        right_element: The name of a GeoDataFrame in `sdata`, or the GeoDataFrame itself
        how: The GeoPandas type of join. By default, left geometries are retained.
        target_coordinate_system: The name of the coordinate system on which the shapes will be transformed. By default, uses the intrinsic coordinate system of the `left_element`.
        **kwargs: Kwargs provided to the [geopandas.sjoin](https://geopandas.org/en/stable/docs/reference/api/geopandas.sjoin.html) function

    Returns:
        The joined `GeoDataFrame`
    """
    if isinstance(left_element, str):
        left_element = sdata[left_element]
    if isinstance(right_element, str):
        right_element = sdata[right_element]

    if target_coordinate_system is None:
        right_element = to_intrinsic(sdata, right_element, left_element)
    else:
        left_element = sdata.transform_element_to_coordinate_system(left_element, target_coordinate_system)
        right_element = sdata.transform_element_to_coordinate_system(right_element, target_coordinate_system)

    return gpd.sjoin(left_element, right_element, how=how, **kwargs)


def _get_cell_id(gdf: gpd.GeoDataFrame, partition: pd.DataFrame, unassigned_value: int | None = 0) -> pd.Series:
    points_gdf = gpd.GeoDataFrame(partition, geometry=gpd.points_from_xy(partition["x"], partition["y"]))
    gdf.index.name = "index_right"  # to reuse the index name later
    spatial_join = points_gdf.sjoin(gdf, how="left")
    spatial_join = spatial_join[~spatial_join.index.duplicated(keep="first")]

    cell_ids = spatial_join["index_right"]

    if unassigned_value is not None:
        cell_ids = (cell_ids.fillna(-1) + 1 + unassigned_value).astype(int)
    else:
        cell_ids = cell_ids.astype("Int64")  # integers with NaNs

    return cell_ids


def assign_transcript_to_cell(
    sdata: SpatialData,
    points_key: str | None = None,
    shapes_key: str | None = None,
    key_added: str = "cell_index",
    unassigned_value: int | None = None,
):
    """Assign each transcript to a cell based on the cell boundaries. It updates the transcript dataframe by adding a new column.

    Args:
        sdata: A `SpatialData` object
        points_key: Key of the spatialdata object containing the transcript dataframe.
        shapes_key: Key of the spatialdata object containing the cell boundaries.
        key_added: Key that will be added to the transcript dataframe containing the cell ID
        unassigned_value: If `None`, transcripts that are not inside any cell will be assigned to NaN. If an integer, this value will be used as the unassigned value.
    """
    df = get_spatial_element(sdata.points, points_key or sdata.attrs.get(SopaAttrs.TRANSCRIPTS))
    geo_df = get_boundaries(sdata, key=shapes_key)

    geo_df = to_intrinsic(sdata, geo_df, df)
    geo_df = geo_df.reset_index()

    get_cell_id = partial(_get_cell_id, geo_df, unassigned_value=unassigned_value)

    if isinstance(df, dd.DataFrame):
        df[key_added] = df.map_partitions(get_cell_id)
    else:
        raise ValueError(f"Invalid dataframe type: {type(df)}")
