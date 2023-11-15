import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import Delaunay

from ._build import _check_has_delaunay
from ._graph import Component

log = logging.getLogger(__name__)


def geometrize_domains(adata: AnnData, domain_key: str) -> gpd.GeoDataFrame:
    _check_has_delaunay(adata)
    data = {"geometry": [], domain_key: []}

    delaunay = Delaunay(adata.obsm["spatial"])
    connectivities = adata.obsp["spatial_connectivities"]
    values = adata.obs[domain_key].values

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
        data[domain_key].append(values[component.first_vertex_index()])

    gdf = gpd.GeoDataFrame(data)

    gdf["length"] = gdf.length
    gdf["area"] = gdf.area
    gdf["roundness"] = 4 * np.pi * gdf["area"] / gdf["length"] ** 2

    return gdf


def domains_geometry_stats(
    adata: AnnData,
    domain_key: str,
    aggregation: str | list[str] = "min",
    key_added_suffix: str = "_distance_to_domain_",
) -> gpd.GeoDataFrame:
    gdf = geometrize_domains(adata, domain_key)

    pairwise_distances: pd.DataFrame = gdf.geometry.apply(lambda g: gdf.distance(g))
    pairwise_distances[domain_key] = gdf[domain_key]

    if isinstance(aggregation, str):
        aggregation = [aggregation]

    for aggr in aggregation:
        df = pairwise_distances.groupby(domain_key).aggregate(aggr).T
        df.columns = [f"{aggr}{key_added_suffix}{c}" for c in df.columns]
        gdf[df.columns] = df

    return gdf.groupby(domain_key)[gdf.columns[2:]].mean()
