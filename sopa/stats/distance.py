import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from tqdm import tqdm

log = logging.getLogger(__name__)


def cells_to_groups(
    adata: AnnData,
    group_key: str,
    key_added_prefix: str | None = "distance_group_",
    ignore_zeros: bool = False,
) -> pd.DataFrame | None:
    distances_to_groups = {}

    if not adata.obs[group_key].dtype.name == "category":
        log.info(f"Converting adata.obs['{group_key}'] to category")
        adata.obs[group_key] = adata.obs[group_key].astype("category")

    for group_id in tqdm(adata.obs[group_key].cat.categories):
        group_nodes = np.where(adata.obs[group_key] == group_id)[0]

        distances = np.full(adata.n_obs, np.nan)

        if not ignore_zeros:
            distances[group_nodes] = 0
            visited = set(group_nodes)
        else:
            visited = set()

        queue = group_nodes
        current_distance = 0

        while len(queue):
            distances[queue] = current_distance

            neighbors = set(adata.obsp["spatial_connectivities"][queue].indices)
            queue = np.array(list(neighbors - visited))
            visited |= neighbors

            current_distance += 1

        distances_to_groups[group_id] = distances

    df_distances = pd.DataFrame(distances_to_groups, index=adata.obs_names)

    if key_added_prefix is None:
        return df_distances
    adata.obsm[f"{key_added_prefix}{group_key}"] = df_distances


def mean_distance(
    adata: AnnData, group_key: str, target_group_key: str | None = None, ignore_zeros: bool = False
) -> pd.DataFrame:
    target_group_key = group_key if target_group_key is None else target_group_key

    df_distances = cells_to_groups(adata, target_group_key, None, ignore_zeros=ignore_zeros)

    if ignore_zeros:
        df_distances.replace(0, np.nan, inplace=True)

    df_distances[group_key] = adata.obs[group_key]
    return df_distances.groupby(group_key, observed=False).mean()
