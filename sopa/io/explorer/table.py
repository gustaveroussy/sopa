from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import zarr
from anndata import AnnData
from scipy.sparse import csr_matrix

from ._constants import FileNames, cell_categories_attrs
from .utils import explorer_file_path

log = logging.getLogger(__name__)


def write_gene_counts(
    path: str, adata: AnnData, layer: str | None = None, is_dir: bool = True
) -> None:
    """Write a `cell_feature_matrix.zarr.zip` file containing the cell-by-gene transcript counts (i.e., from `adata.X`)

    Args:
        path: Path to the Xenium Explorer directory where the cell-by-gene file will be written
        adata: An `AnnData` object. Note that `adata.X` has to be a sparse matrix (and contain the raw counts), else use the `layer` argument.
        layer: If not `None`, `adata.layers[layer]` should be sparse (and contain the raw counts).
        is_dir: If `False`, then `path` is a path to a single file, not to the Xenium Explorer directory.
    """
    path = explorer_file_path(path, FileNames.TABLE, is_dir)

    log.info(f"Writing table with {adata.n_vars} columns")
    counts = adata.X if layer is None else adata.layers[layer]
    counts = csr_matrix(counts.T)

    feature_keys = list(adata.var_names) + ["Total transcripts"]
    feature_ids = feature_keys
    feature_types = ["gene"] * len(adata.var_names) + ["aggregate_gene"]

    ATTRS = {
        "major_version": 3,
        "minor_version": 0,
        "number_cells": adata.n_obs,
        "number_features": adata.n_vars + 1,
        "feature_keys": feature_keys,
        "feature_ids": feature_ids,
        "feature_types": feature_types,
    }

    total_counts = counts.sum(1).A1
    loc = total_counts > 0

    data = np.concatenate([counts.data, total_counts[loc]])
    indices = np.concatenate([counts.indices, np.where(loc)[0]])
    indptr = counts.indptr
    indptr = np.append(indptr, indptr[-1] + loc.sum())

    cell_id = np.ones((adata.n_obs, 2))
    cell_id[:, 0] = np.arange(adata.n_obs)

    with zarr.ZipStore(path, mode="w") as store:
        g = zarr.group(store=store)
        cells_group = g.create_group("cell_features")
        cells_group.attrs.put(ATTRS)

        cells_group.array("cell_id", cell_id, dtype="uint32", chunks=cell_id.shape)
        cells_group.array("data", data, dtype="uint32", chunks=data.shape)
        cells_group.array("indices", indices, dtype="uint32", chunks=indices.shape)
        cells_group.array("indptr", indptr, dtype="uint32", chunks=indptr.shape)


def _write_categorical_column(
    root: zarr.Group, index: int, values: np.ndarray, categories: list[str]
) -> None:
    group = root.create_group(index)
    values_indices = [np.where(values == cat)[0] for cat in categories]
    values_cum_len = np.cumsum([len(indices) for indices in values_indices])

    indices = np.concatenate(values_indices)
    indptr = np.concatenate([[0], values_cum_len[:-1]])

    group.array("indices", indices, dtype="uint32", chunks=(len(indices),))
    group.array("indptr", indptr, dtype="uint32", chunks=(len(indptr),))


def write_cell_categories(path: str, adata: AnnData, is_dir: bool = True) -> None:
    """Write a `analysis.zarr.zip` file containing the cell categories/clusters (i.e., from `adata.obs`)

    Args:
        path: Path to the Xenium Explorer directory where the cell-categories file will be written
        adata: An `AnnData` object
        is_dir: If `False`, then `path` is a path to a single file, not to the Xenium Explorer directory.
    """
    path = explorer_file_path(path, FileNames.CELL_CATEGORIES, is_dir)

    adata.strings_to_categoricals()
    cat_columns = [name for name, cat in adata.obs.dtypes.items() if cat == "category"]

    log.info(f"Writing {len(cat_columns)} cell categories: {', '.join(cat_columns)}")

    ATTRS = cell_categories_attrs()
    ATTRS["number_groupings"] = len(cat_columns)

    with zarr.ZipStore(path, mode="w") as store:
        g = zarr.group(store=store)
        cell_groups = g.create_group("cell_groups")

        for i, name in enumerate(cat_columns):
            if adata.obs[name].isna().any():
                NA = "NA"
                log.warn(f"Column {name} has nan values. They will be displayed as '{NA}'")
                adata.obs[name] = adata.obs[name].cat.add_categories(NA).fillna(NA)

            categories = list(adata.obs[name].cat.categories)
            ATTRS["grouping_names"].append(name)
            ATTRS["group_names"].append(categories)

            _write_categorical_column(cell_groups, i, adata.obs[name], categories)

        cell_groups.attrs.put(ATTRS)


def save_column_csv(path: str, adata: AnnData, key: str):
    """Save one column of the AnnData object as a CSV that can be open interactively in the explorer, under the "cell" panel.

    Args:
        path: Path where to write the CSV that will be open in the Xenium Explorer
        adata: An `AnnData` object
        key: Key of `adata.obs` containing the column to convert
    """
    df = pd.DataFrame({"cell_id": adata.obs_names, "group": adata.obs[key].values})
    df.to_csv(path, index=None)
