import logging

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel, TableModel
from spatialdata.transformations import get_transformation

from ..._constants import SopaKeys
from ..._sdata import get_intrinsic_cs, get_item, get_key
from .. import shapes
from .aggregate import read_all_baysor_patches

log = logging.getLogger(__name__)


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
    n: int | None = None,
    min_area: float = 0,
):
    patches_cells, adatas = read_all_baysor_patches(baysor_temp_dir, min_area, n)
    geo_df, cells_indices, new_ids = resolve(patches_cells, adatas)

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
