import anndata
import numpy as np
import zarr


def add_group(root: zarr.Group, index: int, values: np.ndarray, categories: list[str]):
    group = root.create_group(index)
    values_indices = [np.where(values == cat)[0] for cat in categories]
    values_cum_len = np.cumsum([len(indices) for indices in values_indices])

    indices = np.concatenate(values_indices)
    indptr = np.concatenate([[0], values_cum_len[:-1]])

    group.array("indices", indices, dtype="uint32", chunks=(len(indices),))
    group.array("indptr", indptr, dtype="uint32", chunks=(len(indptr),))


def write_groups(path: str, adata: anndata.AnnData):
    ATTRS = {
        "major_version": 1,
        "minor_version": 0,
        "number_groupings": 1,
        "grouping_names": [],
        "group_names": [],
    }

    categorical_columns = [
        name for name, cat in adata.obs.dtypes.items() if cat == "category"
    ]

    with zarr.ZipStore(path, mode="w") as store:
        g = zarr.group(store=store)
        cell_groups = g.create_group("cell_groups")

        for i, name in enumerate(categorical_columns):
            categories = list(adata.obs[name].cat.categories)
            ATTRS["grouping_names"].append(name)
            ATTRS["group_names"].append(categories)

            add_group(cell_groups, i, adata.obs[name], categories)

        cell_groups.attrs.put(ATTRS)
