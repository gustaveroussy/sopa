import numpy as np
import pandas as pd
import zarr


def add_group(root: zarr.Group, index: int, categories: np.ndarray):
    group = root.create_group(index)
    categories_indices = [
        np.where(categories == cat)[0] for cat in np.unique(categories)
    ]
    categories_cum_len = np.cumsum([len(indices) for indices in categories_indices])

    indices = np.concatenate(categories_indices)
    indptr = np.concatenate([[0], categories_cum_len[:-1]])

    group.array("indices", indices, dtype="uint32", chunks=(len(indices),))
    group.array("indptr", indptr, dtype="uint32", chunks=(len(indptr),))


def write_groups(path: str, df: pd.DataFrame):
    ATTRS = {
        "major_version": 1,
        "minor_version": 0,
        "number_groupings": 1,
        "grouping_names": [],
        "group_names": [],
    }

    with zarr.ZipStore(path, mode="w") as store:
        g = zarr.group(store=store)
        cell_groups = g.create_group("cell_groups")

        i = 0
        for name in df.columns:
            if df[name].dtype == "category":
                ATTRS["grouping_names"].append(name)
                ATTRS["group_names"].append(list(df[name].unique()))

                add_group(cell_groups, i, df[name])
                i += 1

        cell_groups.attrs.put(ATTRS)
