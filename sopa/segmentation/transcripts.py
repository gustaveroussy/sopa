from __future__ import annotations

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
from tqdm import tqdm

from .._constants import SopaKeys
from .._sdata import get_element, get_key, save_shapes, save_table
from . import aggregate, shapes

log = logging.getLogger(__name__)


def resolve(
    sdata: SpatialData,
    temp_dir: str,
    gene_column: str,
    patches_dirs: list[str] | None = None,
    min_area: float = 0,
    shapes_key: str = SopaKeys.BAYSOR_BOUNDARIES,
):
    """Concatenate all the per-patch segmentation runs and resolve the conflicts. Resulting cells boundaries are saved in the `SpatialData` object.

    Args:
        sdata: A `SpatialData` object
        temp_dir: Temporary directory used to store all the patches subdirectories (one subdirectory for one patch and for one segmentation run)
        gene_column: Column of the transcript dataframe containing the genes names
        patches_dirs: Optional list of subdirectories inside `temp_dir` to be read. By default, read all.
        min_area: Minimum area (in microns^2) for a cell to be kept
    """
    if min_area > 0:
        log.info(f"Cells whose area is less than {min_area} microns^2 will be removed")

    patches_cells, adatas = _read_all_segmented_patches(temp_dir, min_area, patches_dirs)
    geo_df, cells_indices, new_ids = _resolve_patches(patches_cells, adatas)

    image_key = get_key(sdata, "images")
    points = get_element(sdata, "points")
    transformations = get_transformation(points, get_all=True).copy()

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
    table.obs[SopaKeys.REGION_KEY] = pd.Series(shapes_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.SLIDE_KEY] = pd.Series(image_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.INSTANCE_KEY] = geo_df.index

    table = TableModel.parse(
        table,
        region_key=SopaKeys.REGION_KEY,
        region=shapes_key,
        instance_key=SopaKeys.INSTANCE_KEY,
    )

    sdata.shapes[shapes_key] = geo_df
    save_shapes(sdata, shapes_key, overwrite=True)

    sdata.tables[SopaKeys.TABLE] = table
    save_table(sdata, SopaKeys.TABLE)

    log.info(
        f"Added sdata.tables['{SopaKeys.TABLE}'], and {len(geo_df)} cell boundaries to sdata['{shapes_key}']"
    )


def _read_one_segmented_patch(
    directory: str, min_area: float = 0, min_vertices: int = 4
) -> tuple[list[Polygon], AnnData]:
    directory: Path = Path(directory)

    loom_file = directory / "segmentation_counts.loom"
    if loom_file.exists():
        adata = anndata.read_loom(
            directory / "segmentation_counts.loom", obs_names="Name", var_names="Name"
        )
    else:
        adata = anndata.read_h5ad(directory / "segmentation_counts.h5ad")

    adata.obs.rename(columns={"area": SopaKeys.ORIGINAL_AREA_OBS}, inplace=True)

    cells_num = pd.Series(adata.obs["CellID"].astype(int), index=adata.obs_names)
    del adata.obs["CellID"]

    with open(directory / "segmentation_polygons.json") as f:
        polygons_dict = json.load(f)
        polygons_dict = {c["cell"]: c for c in polygons_dict["geometries"]}

    cells_num = cells_num[
        cells_num.map(lambda num: len(polygons_dict[num]["coordinates"][0]) >= min_vertices)
    ]

    gdf = gpd.GeoDataFrame(
        index=cells_num.index, geometry=[shape(polygons_dict[cell_num]) for cell_num in cells_num]
    )

    gdf.geometry = gdf.geometry.map(lambda cell: shapes._ensure_polygon(cell))
    gdf = gdf[~gdf.geometry.isna()]

    ratio_filtered = (gdf.area <= min_area).mean()
    if ratio_filtered > 0.2:
        log.warn(f"{ratio_filtered:.2%} of cells will be filtered due to {min_area=}")

    gdf = gdf[gdf.area > min_area]

    return gdf.geometry.values, adata[gdf.index].copy()


def _read_all_segmented_patches(
    temp_dir: str,
    min_area: float = 0,
    patches_dirs: list[str] | None = None,
) -> tuple[list[list[Polygon]], list[AnnData]]:
    if patches_dirs is None or not len(patches_dirs):
        patches_dirs = [subdir for subdir in Path(temp_dir).iterdir() if subdir.is_dir()]

    outs = [
        _read_one_segmented_patch(path, min_area)
        for path in tqdm(patches_dirs, desc="Reading transcript-segmentation outputs")
    ]

    patches_cells, adatas = zip(*outs)

    return patches_cells, adatas


def _resolve_patches(
    patches_cells: list[list[Polygon]], adatas: list[AnnData]
) -> tuple[gpd.GeoDataFrame, np.ndarray, np.ndarray]:
    """Resolve the segmentation conflits on the patches overlaps.

    Args:
        patches_cells: List of polygons segmented on each patch
        adatas: List of AnnData objects corresponding to each patch

    Returns:
        The new GeoDataFrame, the new cells indices (-1 for merged cells), and the ids of the merged cells.
    """
    patch_ids = [adata.obs_names for adata in adatas]

    patch_indices = np.arange(len(patches_cells)).repeat([len(cells) for cells in patches_cells])
    cells = [cell for cells in patches_cells for cell in cells]
    segmentation_ids = np.array([cell_id for ids in patch_ids for cell_id in ids])

    cells_resolved, cells_indices = shapes.solve_conflicts(
        cells, patch_indices=patch_indices, return_indices=True
    )

    existing_ids = segmentation_ids[cells_indices[cells_indices >= 0]]
    new_ids = np.char.add("merged_cell_", np.arange((cells_indices == -1).sum()).astype(str))
    index = np.concatenate([existing_ids, new_ids])

    return (
        gpd.GeoDataFrame({"geometry": cells_resolved}, index=index),
        cells_indices,
        new_ids,
    )


def copy_segmentation_config(path: Path, config: dict, config_path: str | None):
    """Copy the segmentation config to a file.

    Args:
        path: Where the config will be saved
        config: Dictionnary config
        config_path: Already existing config file, will be copied if provided
    """
    if config_path is not None:
        import shutil

        shutil.copy(config_path, path)
        return

    if path.suffix == ".json":
        with open(path, "w") as f:
            json.dump(config, f)
            return

    if path.suffix == ".toml":
        try:
            import toml
        except ImportError:
            raise ImportError(
                "To use baysor, you need its corresponding sopa extra: `pip install 'sopa[baysor]'` (normal mode) or `pip install -e '.[baysor]'` (if using snakemake).\
                \nAlso, make sure to install the baysor executable (https://github.com/kharchenkolab/Baysor)."
            )

        with open(path, "w") as f:
            toml.dump(config, f)
            return

    raise ValueError(f"Config file must be either a .json or a .toml file. Found: {path.suffix}")
