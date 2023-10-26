import numpy as np
import pandas as pd
from anndata import AnnData


def cells_to_domains(
    adata: AnnData, domain_key: str, key_added_prefix: str | None = "distance_domain_"
) -> pd.DataFrame | None:
    distances_to_domains = {}

    for domain_id in adata.obs[domain_key].unique():
        domain_nodes = np.where(adata.obs[domain_key] == domain_id)[0]

        distances = np.full(adata.n_obs, np.nan)
        current_distance = 0
        distances[domain_nodes] = current_distance

        visited = set(domain_nodes)
        queue = domain_nodes

        while len(queue):
            distances[queue] = current_distance

            neighbors = set(adata.obsp["spatial_connectivities"][queue].indices)
            queue = np.array(list(neighbors - visited))
            visited |= neighbors

            current_distance += 1

        if key_added_prefix is None:
            distances_to_domains[domain_id] = distances
        else:
            adata.obs[f"{key_added_prefix}{domain_id}"] = distances

    if key_added_prefix is None:
        return pd.DataFrame(distances_to_domains, index=adata.obs_names)


def mean_distance_group_to_group():
    ...


def mean_distance_domain_to_domain():
    ...


def mean_distance_group_to_domains(adata: AnnData, domain_key: str, group_key: str) -> pd.DataFrame:
    series = []
    for domain_id in adata.obs[domain_key].unique():
        key = f"distance_to_niche_{domain_id}"
        obs = adata.obs
        obs = obs[obs[key] > 0]
        s = obs.groupby(group_key)[key].mean()
        series.append(s)

    return pd.concat(
        series,
        axis=1,
        keys=[f"Mean distance to niche: {niche}" for niche in adata.obs[domain_key].unique()],
    )
