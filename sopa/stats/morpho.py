import logging
from collections import defaultdict

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse.csgraph import connected_components
from tqdm import tqdm

log = logging.getLogger(__name__)


def triangle_area(adata: AnnData, i: int, j: int, k: int) -> float:
    coordinates = adata.obsm["spatial"][[i, j, k]]
    dy = np.roll(coordinates[:, 1], -1) - np.roll(coordinates[:, 1], 1)
    return np.abs((coordinates[:, 0] * dy).sum()) / 2


def _component_morphology(adata_component: AnnData) -> list[float]:
    dist = adata_component.obsp["spatial_distances"]

    to_visit = {0}
    visited = set()

    triangle_counter = defaultdict(int)  # pairs of edge indices
    triangles = set()  # tuple of 3 indices

    while to_visit:
        index = to_visit.pop()
        visited.add(index)

        nghs = dist[index].indices
        ngh_i, ngh_j = dist[nghs][:, nghs].nonzero()
        ngh_i, ngh_j = nghs[ngh_i], nghs[ngh_j]

        for i, j in zip(ngh_i, ngh_j):
            sorted_indices = tuple(sorted((index, i, j)))
            triangles.add(sorted_indices)
            triangle_counter[(sorted_indices[0], sorted_indices[1])] += 1
            triangle_counter[(sorted_indices[0], sorted_indices[2])] += 1
            triangle_counter[(sorted_indices[1], sorted_indices[2])] += 1

        to_visit |= set(nghs) - visited

    perimeter = sum(dist[i, j] for (i, j), count in triangle_counter.items() if count == 6)
    area = sum(triangle_area(adata_component, *indices) for indices in triangles)
    roundness = 4 * np.pi * area / perimeter**2

    return {"perimeter": perimeter, "area": area, "roundness": roundness}


def morphology(
    adata: AnnData, domain_key: str, domain_id: str, min_component_size: int = 20
) -> pd.DataFrame:
    mask = adata.obs[domain_key] == domain_id
    subgraph = adata.obsp["spatial_connectivities"][mask][:, mask]
    n_components, component_labels = connected_components(subgraph, directed=False)

    component_counts = 0
    geo_stats = []

    for label in tqdm(range(n_components)):
        component_indices = np.where(mask)[0][np.where(component_labels == label)[0]]
        if len(component_indices) >= min_component_size:
            component_counts += 1

            adata_ = adata[component_indices].copy()
            geo_stats.append(_component_morphology(adata_))

    if geo_stats:
        res = pd.concat([pd.Series(d) for d in geo_stats], axis=1).mean(1)
        res["count"] = component_counts
    else:
        res = pd.Series([np.nan, np.nan, np.nan, np.nan])
        res.index = ["perimeter", "area", "roundness", "count"]

    res.name = domain_id
    return res


def morphologies(adata: AnnData, domain_key: str, min_component_size: int = 20):
    domain_ids = adata.obs[domain_key].unique()
    log.info(f"Computing morphologies of {len(domain_ids)} domains")
    res = pd.concat(
        [
            morphology(adata, domain_key, domain_id, min_component_size=min_component_size)
            for domain_id in domain_ids
        ],
        axis=1,
    ).T
    res.index.name = domain_key
    return res


def interface_ratio(adata: AnnData, domain_key: str):
    ...
