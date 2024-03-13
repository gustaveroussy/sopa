from __future__ import annotations

import geopandas as gpd
from spatialdata import SpatialData

from .._sdata import get_intrinsic_cs


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
        target_coordinate_system = get_intrinsic_cs(sdata, left_element)

    left_element = sdata.transform_element_to_coordinate_system(
        left_element, target_coordinate_system
    )
    right_element = sdata.transform_element_to_coordinate_system(
        right_element, target_coordinate_system
    )

    return gpd.sjoin(left_element, right_element, how=how, **kwargs)
