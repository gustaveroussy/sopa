from __future__ import annotations

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import Delaunay
from shapely.geometry import MultiPolygon, Polygon
from spatialdata import SpatialData

from .._constants import SopaKeys
from ._build import _check_has_delaunay
from ._graph import Component

log = logging.getLogger(__name__)


def geometrize_niches(
    adata: AnnData | SpatialData,
    niche_key: str,
    buffer: int | str = "auto",
    perc_area_th: float = 0.05,
) -> gpd.GeoDataFrame:
    """Converts the niches to shapely polygons, and put into a `GeoDataFrame`. Note that each niche can appear multiple times, as they can be separated by other niches ; in this case, we call them different "components" of the same niche ID.

    Plot components:
        You can show niches components with GeoPandas
        ```py
        gdf = geometrize_niches(adata, niche_key)
        gdf.plot(column=niche_key)
        ```

    Args:
        adata: An `AnnData` object, or a `SpatialData object`
        niche_key: Key of `adata.obs` containing the niches
        buffer: Expansion radius applied on components. By default, `3 * mean_distance_neighbors`
        perc_area_th: For each niche, components whose area is less than `perc_area_th * max_component_area` will be removed

    Returns:
        A `GeoDataFrame` with geometries for each niche component. We also compute the area/perimeter/roundness of each component.
    """
    if isinstance(adata, SpatialData):
        adata = adata.tables[SopaKeys.TABLE]

    _check_has_delaunay(adata)
    data = {"geometry": [], niche_key: []}

    delaunay = Delaunay(adata.obsm["spatial"])
    connectivities = adata.obsp["spatial_connectivities"]
    values = adata.obs[niche_key].values

    keep = (
        (connectivities[delaunay.simplices[:, 0], delaunay.simplices[:, 1]].A1 == 1)
        & (connectivities[delaunay.simplices[:, 0], delaunay.simplices[:, 2]].A1 == 1)
        & (connectivities[delaunay.simplices[:, 1], delaunay.simplices[:, 2]].A1 == 1)
        & (values[delaunay.simplices[:, 0]] == values[delaunay.simplices[:, 1]])
        & (values[delaunay.simplices[:, 0]] == values[delaunay.simplices[:, 2]])
        & (values[delaunay.simplices[:, 0]] == values[delaunay.simplices[:, 2]])
    )  # Keep simplices that are in the original Delaunay graph, and which are not in between different value categories

    neighbors = np.where(np.isin(delaunay.neighbors, np.where(~keep)[0]), -1, delaunay.neighbors)

    simplices_to_visit = set(np.where(keep)[0])

    while simplices_to_visit:
        component = Component(adata, delaunay, neighbors)
        component.visit(simplices_to_visit)

        data["geometry"].append(component.polygon)
        data[niche_key].append(values[component.first_vertex_index()])

    gdf = gpd.GeoDataFrame(data)

    if buffer is not None and buffer != 0:
        gdf = _clean_components(adata, gdf, niche_key, buffer)

    gdf[SopaKeys.GEOMETRY_LENGTH] = gdf.length
    gdf[SopaKeys.GEOMETRY_AREA] = gdf.area
    gdf[SopaKeys.GEOMETRY_ROUNDNESS] = (
        4 * np.pi * gdf[SopaKeys.GEOMETRY_AREA] / gdf[SopaKeys.GEOMETRY_LENGTH] ** 2
    )

    # Remove minor components (compared to the largest component of its corresponding niche)
    gdf = gdf[gdf.area >= gdf[niche_key].map(gdf.groupby(niche_key).area.max() * perc_area_th)]

    return gdf


def _clean_components(
    adata: AnnData, gdf: gpd.GeoDataFrame, niche_key: str, buffer: int | str
) -> gpd.GeoDataFrame:
    data = {"geometry": [], niche_key: []}

    if buffer == "auto":
        buffer = 3 * adata.obsp["spatial_distances"].data.mean()

    for niche, group_gdf in gdf.groupby(niche_key):
        multi_polygon = MultiPolygon(group_gdf.geometry.values).buffer(buffer).buffer(-buffer)

        if isinstance(multi_polygon, Polygon):
            data["geometry"].append(multi_polygon)
            data[niche_key].append(niche)
            continue

        for geometry in multi_polygon.geoms:
            data["geometry"].append(geometry)
            data[niche_key].append(niche)

    return gpd.GeoDataFrame(data)


def niches_geometry_stats(
    adata: AnnData | SpatialData,
    niche_key: str,
    aggregation: str | list[str] = "min",
    key_added_suffix: str = "_distance_to_niche_",
    **geometrize_niches_kwargs: str,
) -> gpd.GeoDataFrame:
    """Computes statistics over niches geometries

    Details:
        - `n_components`: Number of connected component of a niche (a component is a group of neighbor cells with the same niche attribute)
        - `length`: Mean distance of the exterior/boundary of the components of a niche
        - `area`: Mean area of the components of a niche
        - `roundness`: Float value between 0 and 1. The higher the value, the closer to a circle. Computed via `4 * pi * area / length**2`
        - `mean_distance_to_niche_X`: mean distance to the niche (between the two closest points of the niches)

    Args:
        adata: An `AnnData` object, or a `SpatialData object`
        niche_key: Key of `adata.obs` containing the niches
        aggregation: Aggregation mode. Either one string such as `"min"`, or a list such as `["mean", "min"]`.
        key_added_suffix: Suffix added in the DataFrame columns. Defaults to "_distance_to_niche_".
        geometrize_niches_kwargs: Kwargs to the `sopa.spatial.geometrize_niches` function

    Returns:
        A `DataFrame` of shape `n_niches * n_statistics`
    """
    if isinstance(adata, SpatialData):
        adata = adata.tables[SopaKeys.TABLE]

    gdf = geometrize_niches(adata, niche_key, **geometrize_niches_kwargs)
    value_counts = gdf[niche_key].value_counts()

    assert len(gdf), "No niche geometry found, stats can't be computed"

    log.info(f"Computing pairwise distances between {len(gdf)} components")
    pairwise_distances: pd.DataFrame = gdf.geometry.apply(lambda g: gdf.distance(g))
    pairwise_distances[niche_key] = gdf[niche_key]

    if isinstance(aggregation, str):
        aggregation = [aggregation]

    for aggr in aggregation:
        df = pairwise_distances.groupby(niche_key).aggregate(aggr).T
        df.columns = [f"{aggr}{key_added_suffix}{c}" for c in df.columns]
        gdf[df.columns] = df

    df_stats: pd.DataFrame = gdf.groupby(niche_key)[gdf.columns[2:]].mean()
    df_stats.insert(0, SopaKeys.GEOMETRY_COUNT, value_counts)

    return df_stats
