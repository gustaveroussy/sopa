import json
import logging
from pathlib import Path

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from shapely.geometry import Polygon, shape
from spatialdata import SpatialData
from spatialdata.models import ShapesModel, TableModel
from spatialdata.transformations import get_transformation

from ..._constants import SopaKeys
from ..._sdata import get_intrinsic_cs, get_item, get_key
from .. import shapes

log = logging.getLogger(__name__)


def read_baysor(
    directory: str, min_area: float = 0, expand_radius: float = 0, min_vertices: int = 4
) -> tuple[list[Polygon], AnnData]:
    directory = Path(directory)

    adata = anndata.read_loom(
        directory / "segmentation_counts.loom", obs_names="Name", var_names="Name"
    )
    adata = adata[adata.obs.area > min_area].copy()

    cells_num = pd.Series(adata.obs_names.str.split("-").str[-1].astype(int), index=adata.obs_names)

    with open(directory / "segmentation_polygons.json") as f:
        polygons_dict = json.load(f)
        polygons_dict = {c["cell"]: c for c in polygons_dict["geometries"]}

    cells_num = cells_num[
        cells_num.map(lambda num: len(polygons_dict[num]["coordinates"][0]) >= min_vertices)
    ]

    cells = [shape(polygons_dict[cell_num]).buffer(expand_radius) for cell_num in cells_num]

    for cell in cells:
        if not isinstance(cell, Polygon):
            print(directory)

    return cells, adata[cells_num.index].copy()


def read_all_baysor_patches(
    baysor_temp_dir: str,
    min_area: float = 0,
    expand_radius: float = 0,
    patches_dirs: list[str] | None = None,
) -> tuple[list[list[Polygon]], list[AnnData]]:
    if patches_dirs is None:
        baysor_temp_dir = Path(baysor_temp_dir)
        outs = [
            read_baysor(directory, min_area, expand_radius)
            for directory in baysor_temp_dir.iterdir()
        ]
    else:
        outs = [read_baysor(path, min_area, expand_radius) for path in patches_dirs]

    patches_cells, adatas = zip(*outs)

    return patches_cells, adatas


def resolve_patches(
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


def resolve(
    sdata: SpatialData,
    baysor_temp_dir: str,
    gene_column: str,
    patches_dirs: list[str] | None = None,
    min_area: float = 0,
    expand_radius: float = 0,
):
    patches_cells, adatas = read_all_baysor_patches(
        baysor_temp_dir, min_area, expand_radius, patches_dirs
    )
    geo_df, cells_indices, new_ids = resolve_patches(patches_cells, adatas)

    image_key = get_key(sdata, "images")
    points_key, points = get_item(sdata, "points")
    transformations = get_transformation(points, get_all=True)

    geo_df = ShapesModel.parse(geo_df, transformations=transformations)

    table_conflicts = []
    if len(new_ids):
        new_cells = geo_df.geometry[cells_indices == -1]
        geo_df_new = gpd.GeoDataFrame({"geometry": new_cells})
        geo_df_new = ShapesModel.parse(geo_df_new, transformations=transformations)

        table_conflicts = sdata.aggregate(
            values=points_key,
            by=geo_df_new,
            value_key=gene_column,
            agg_func="count",
            target_coordinate_system=get_intrinsic_cs(sdata, points_key),
        ).table
        table_conflicts.obs_names = new_ids
        table_conflicts = [table_conflicts]

    valid_ids = set(list(geo_df.index))
    table = anndata.concat(
        [adata[list(valid_ids & set(list(adata.obs_names)))] for adata in adatas] + table_conflicts,
        join="outer",
    )
    table.obs.dropna(axis="columns", inplace=True)

    geo_df = geo_df.loc[table.obs_names]

    table.obsm["spatial"] = np.array([[centroid.x, centroid.y] for centroid in geo_df.centroid])
    table.obs[SopaKeys.REGION_KEY] = pd.Series(
        SopaKeys.BAYSOR_BOUNDARIES, index=table.obs_names, dtype="category"
    )
    table.obs[SopaKeys.SLIDE_KEY] = pd.Series(image_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.INSTANCE_KEY] = geo_df.index

    table = TableModel.parse(
        table,
        region_key=SopaKeys.REGION_KEY,
        region=SopaKeys.BAYSOR_BOUNDARIES,
        instance_key=SopaKeys.INSTANCE_KEY,
    )

    sdata.add_shapes(SopaKeys.BAYSOR_BOUNDARIES, geo_df, overwrite=True)

    if sdata.table is not None:
        del sdata.table

    sdata.table = table

    log.info(
        f"Added sdata.table, and {len(geo_df)} cell boundaries to sdata['{SopaKeys.BAYSOR_BOUNDARIES}']"
    )
