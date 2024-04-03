from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialdata import SpatialData
from tqdm import tqdm

from .._constants import SopaKeys
from ._build import _check_has_delaunay

log = logging.getLogger(__name__)


def cells_to_groups(
    adata: AnnData,
    group_key: str,
    key_added_prefix: str | None = None,
    ignore_zeros: bool = False,
) -> pd.DataFrame | None:
    """Compute the hop-distance between each cell and a cell category/group.

    Args:
        adata: An `AnnData` object, or a `SpatialData object`
        group_key: Key of `adata.obs` containing the groups
        key_added_prefix: Prefix to the key added in `adata.obsm`. If `None`, will return the `DataFrame` instead of saving it.
        ignore_zeros: If `True`, a cell distance to its own group is 0.

    Returns:
        A `Dataframe` of shape `n_obs * n_groups`, or `None` if `key_added_prefix` was provided (in this case, the dataframe is saved in `"{key_added_prefix}{group_key}"`)
    """
    if isinstance(adata, SpatialData):
        adata = adata.tables[SopaKeys.TABLE]

    _check_has_delaunay(adata)

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
    """Mean distance between two groups (typically, between cell-types, or between cell-types and domains)

    Note:
        The distance is a number of hops, i.e. a distance of 10 between a pDC and a T cell means that there are 10 cells on the closest path from one to the other cell.

    Args:
        adata: An `AnnData` object, or a `SpatialData object`
        group_key: Key of `adata.obs` containing the groups
        target_group_key: Key of `adata.obs` containing the target groups (by default, uses `group_key`)
        ignore_zeros: If `True`, a cell distance to its own group is 0.

    Returns:
        `DataFrame` of shape `n_groups * n_groups_target` of mean hop-distances
    """
    if isinstance(adata, SpatialData):
        adata = adata.tables[SopaKeys.TABLE]

    target_group_key = group_key if target_group_key is None else target_group_key

    df_distances = cells_to_groups(adata, target_group_key, None, ignore_zeros=ignore_zeros)

    if ignore_zeros:
        df_distances.replace(0, np.nan, inplace=True)

    df_distances[group_key] = adata.obs[group_key]
    df_distances = df_distances.groupby(group_key, observed=False).mean()
    df_distances.columns.name = target_group_key

    return df_distances


def prepare_network(
    adata: AnnData,
    cell_type_key: str,
    niche_key: str,
    clip_weight: float = 3,
    node_colors: tuple[str] = ("#5c7dc4", "#f05541"),
    node_sizes: float = (1.3, 5),
) -> tuple[pd.DataFrame, dict, dict, dict]:
    """Create a dataframe representing weights between cell-types and/or niches.
    This can be later use to plot a cell-type/niche represention of a whole slide
    using the netgraph library.

    Args:
        adata: An `AnnData` object
        cell_type_key: Key of `adata.obs` containing the cell types
        niche_key: Key of `adata.obs` containing the niches
        clip_weight: Maximum weight
        node_colors: Tuple of (cell-type color, niche color)
        node_sizes: Tuple of (cell-type size, niche size)

    Returns:
        A DataFrame of weights between cell-types and/or niches, and three dict for netgraph display
    """
    node_color, node_size, node_shape = {}, {}, {}

    log.info("Computing all distances for the 4 pairs of categories")
    weights = mean_distance(adata, cell_type_key)
    top_right = mean_distance(adata, cell_type_key, niche_key)
    bottom_left = mean_distance(adata, niche_key, cell_type_key)
    bottom_right = mean_distance(adata, niche_key, niche_key)

    for pop in weights.index:
        node_color[pop] = node_colors[0]
        node_size[pop] = node_sizes[0]
        node_shape[pop] = "o"

    for niche in bottom_right.index:
        node_color[niche] = node_colors[1]
        node_size[niche] = node_sizes[1]
        node_shape[niche] = "h"

    # assemble dataframe per-block
    bottom_left[bottom_right.columns] = bottom_right
    weights[top_right.columns] = top_right
    weights = pd.concat([weights, bottom_left], axis=0).copy()

    # convert distances to symmetric weights
    weights = 1 / weights
    np.fill_diagonal(weights.values, 0)
    weights = weights.clip(0, clip_weight)
    weights = (weights.T + weights) / 2

    return weights, node_color, node_size, node_shape
