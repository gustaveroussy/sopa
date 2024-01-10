import json
import logging
from pathlib import Path

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from shapely.geometry import Polygon, shape
from shapely.validation import make_valid
from spatialdata import SpatialData
from spatialdata.models import ShapesModel, TableModel
from spatialdata.transformations import get_transformation
from tqdm import tqdm

from ..._constants import SopaKeys
from ..._sdata import get_element, get_key
from .. import aggregate, shapes

log = logging.getLogger(__name__)


def read_baysor(
    directory: str, min_area: float = 0, min_vertices: int = 4
) -> tuple[list[Polygon], AnnData]:
    directory = Path(directory)

    adata = anndata.read_loom(
        directory / "segmentation_counts.loom", obs_names="Name", var_names="Name"
    )
    adata = adata[adata.obs.area > min_area]

    cells_num = pd.Series(adata.obs_names.str.split("-").str[-1].astype(int), index=adata.obs_names)

    with open(directory / "segmentation_polygons.json") as f:
        polygons_dict = json.load(f)
        polygons_dict = {c["cell"]: c for c in polygons_dict["geometries"]}

    cells_num = cells_num[
        cells_num.map(lambda num: len(polygons_dict[num]["coordinates"][0]) >= min_vertices)
    ]

    gdf = gpd.GeoDataFrame(
        index=cells_num.index, geometry=[shape(polygons_dict[cell_num]) for cell_num in cells_num]
    )

    gdf.geometry = gdf.geometry.map(lambda cell: shapes._ensure_polygon(make_valid(cell)))
    gdf = gdf[~gdf.geometry.isna()]

    return gdf.geometry.values, adata[gdf.index].copy()


def read_all_baysor_patches(
    baysor_temp_dir: str,
    min_area: float = 0,
    patches_dirs: list[str] = None,
) -> tuple[list[list[Polygon]], list[AnnData]]:
    if patches_dirs is None or not len(patches_dirs):
        patches_dirs = [subdir for subdir in Path(baysor_temp_dir).iterdir() if subdir.is_dir()]

    outs = [
        read_baysor(path, min_area) for path in tqdm(patches_dirs, desc="Reading baysor outputs")
    ]

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
    patches_dirs: list[str],
    min_area: float = 0,
):
    """Concatenate all the per-patch Baysor run and resolve the conflicts. Resulting cells boundaries are saved in the `SpatialData` object.

    Args:
        sdata: A `SpatialData` object
        baysor_temp_dir: Temporary directory used to store all the baysor subdirectories (one subdirectory for one patch and for one baysor run)
        gene_column: Column of the transcript dataframe containing the genes names
        patches_dirs: Optional list of subdirectories inside `baysor_temp_dir` to be read. By default, read all.
        min_area: Minimum area (in microns^2) for a cell to be kept
    """
    patches_cells, adatas = read_all_baysor_patches(baysor_temp_dir, min_area, patches_dirs)
    geo_df, cells_indices, new_ids = resolve_patches(patches_cells, adatas)

    image_key = get_key(sdata, "images")
    points = get_element(sdata, "points")
    transformations = get_transformation(points, get_all=True)

    geo_df = ShapesModel.parse(geo_df, transformations=transformations)

    table_conflicts = []
    if len(new_ids):
        new_cells = geo_df.geometry[cells_indices == -1]
        geo_df_new = gpd.GeoDataFrame({"geometry": new_cells})
        geo_df_new = ShapesModel.parse(geo_df_new, transformations=transformations)

        log.info("Aggregating transcripts on merged cells")
        table_conflicts = aggregate.count_transcripts(sdata, gene_column, geo_df=geo_df_new)
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
