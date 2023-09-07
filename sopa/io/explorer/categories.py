import numpy as np
import pandas as pd
import zarr


def add_group(root: zarr.Group, index: int, values: np.ndarray, categories: list[str]):
    group = root.create_group(index)
    values_indices = [np.where(values == cat)[0] for cat in categories]
    values_cum_len = np.cumsum([len(indices) for indices in values_indices])

    indices = np.concatenate(values_indices)
    indptr = np.concatenate([[0], values_cum_len[:-1]])

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
                categories = list(df[name].cat.categories)
                ATTRS["grouping_names"].append(name)
                ATTRS["group_names"].append(categories)

                add_group(cell_groups, i, df[name], categories)
                i += 1

        cell_groups.attrs.put(ATTRS)
