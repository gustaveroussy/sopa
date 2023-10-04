import json
from pathlib import Path

import anndata
import pandas as pd
from anndata import AnnData
from shapely.geometry import Polygon, shape


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
    baysor_temp_dir: str, min_area: float = 0, n: int | None = None
) -> tuple[list[list[Polygon]], list[AnnData]]:
    baysor_temp_dir = Path(baysor_temp_dir)

    if n is None:
        outs = [read_baysor(directory, min_area) for directory in baysor_temp_dir.iterdir()]
    else:
        outs = [read_baysor(baysor_temp_dir / str(i), min_area) for i in range(n)]

    patches_cells, adatas = zip(*outs)

    return patches_cells, adatas
