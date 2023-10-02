import json
from pathlib import Path

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from shapely.geometry import Polygon, shape

from .. import shapes


def read_baysor(
    directory: str, min_area: float = 0, min_vertices: int = 4
) -> tuple[list[Polygon], AnnData]:
    directory = Path(directory)

    adata = anndata.read_loom(
        directory / "segmentation_counts.loom", obs_names="Name", var_names="Name"
    )
    adata = adata[adata.obs.area >= min_area].copy()

    cells_num = pd.Series(adata.obs_names.str.split("-").str[-1].astype(int), index=adata.obs_names)

    with open(directory / "segmentation_polygons.json") as f:
        polygons_dict = json.load(f)
        polygons_dict = {c["cell"]: c for c in polygons_dict["geometries"]}

    cells_num = cells_num[
        cells_num.map(lambda num: len(polygons_dict[num]["coordinates"][0]) >= min_vertices)
    ]

    cells = [shape(polygons_dict[cell_num]) for cell_num in cells_num]

    return cells, adata[cells_num.index].copy()


def read_all_baysor_patches(
    baysor_dir: str, min_area: float = 0, n: int = None
) -> tuple[list[list[Polygon]], list[AnnData]]:
    baysor_dir = Path(baysor_dir)

    if n is None:
        outs = [read_baysor(directory, min_area) for directory in baysor_dir.iterdir()]
    else:
        outs = [read_baysor(baysor_dir / str(i), min_area) for i in range(n)]

    patches_cells, adatas = zip(*outs)

    return patches_cells, adatas


def resolve(
    patches_cells: list[list[Polygon]], adatas: list[AnnData]
) -> tuple[gpd.GeoDataFrame, np.ndarray, np.ndarray]:
    patch_ids = [adata.obs_names for adata in adatas]

    patch_indices = np.arange(len(patches_cells)).repeat([len(cells) for cells in patches_cells])
    cells = [cell for cells in patches_cells for cell in cells]
    baysor_ids = np.array([cell_id for ids in patch_ids for cell_id in ids])

    cells_resolved, cells_indices = shapes.solve_conflicts(
        cells, patch_indices=patch_indices, return_indices=True
    )

    existing_ids = baysor_ids[cells_indices[cells_indices >= 0]]
    new_ids = np.char.add("merged_cell_", np.arange((cells_indices == -1).sum()).astype(str))
    index = np.concatenate([existing_ids, new_ids])

    return (
        gpd.GeoDataFrame({"geometry": cells_resolved}, index=index),
        cells_indices,
        new_ids,
    )
