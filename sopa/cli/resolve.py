import typer

app_resolve = typer.Typer()
option = typer.Option()


@app_resolve.command()
def cellpose(
    sdata_path: str,
    patch_dir: str = option,
    expand_radius: float = typer.Option(default=0),
):
    from pathlib import Path

    import spatialdata
    import zarr
    from shapely.geometry import Polygon
    from tqdm import tqdm

    from sopa.segmentation import shapes
    from sopa.segmentation.update import update
    from sopa.utils.utils import _get_spatial_image

    sdata = spatialdata.read_zarr(sdata_path)

    image_key, _ = _get_spatial_image(sdata)

    polygons = []

    files = [f for f in Path(patch_dir).iterdir() if f.suffix == ".zip"]
    for file in tqdm(files):
        z = zarr.open(file, mode="r")
        for _, coords_zarr in z.arrays():
            polygons.append(Polygon(coords_zarr[:]))
    polygons = shapes.solve_conflicts(polygons)
    polygons = shapes.expand(polygons, expand_radius)

    update(sdata, polygons, image_key)


@app_resolve.command()
def baysor(
    sdata_path: str,
    baysor_dir: str = option,
    gene_column: str = option,
    min_area: float = 0,
    n: int = None,
):
    import anndata
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    import spatialdata
    from spatialdata.models import ShapesModel, TableModel

    from sopa.segmentation.baysor.aggregate import read_all_baysor_patches, resolve
    from sopa.utils.utils import get_intrinsic_cs

    sdata = spatialdata.read_zarr(sdata_path)

    patch_polygons, adatas = read_all_baysor_patches(baysor_dir, min_area, n)
    geo_df, polys_indices, new_ids = resolve(patch_polygons, adatas)
    geo_df = ShapesModel.parse(geo_df, transformations=sdata["transcripts"].attrs["transform"])

    table_conflicts = []
    if len(new_ids):
        new_polys = geo_df.geometry[polys_indices == -1]
        geo_df_new = gpd.GeoDataFrame({"geometry": new_polys})
        geo_df_new = ShapesModel.parse(
            geo_df_new, transformations=sdata["transcripts"].attrs["transform"]
        )

        table_conflicts = sdata.aggregate(
            values="transcripts",
            by=geo_df_new,
            value_key=gene_column,
            agg_func="count",
            target_coordinate_system=get_intrinsic_cs(sdata, "transcripts"),
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
    table.obs["region"] = pd.Series("polygons", index=table.obs_names, dtype="category")
    # table.obs["slide"] = pd.Series(image_key, index=table.obs_names, dtype="category")
    # table.obs["dataset_id"] = pd.Series(image_key, index=table.obs_names, dtype="category")
    table.obs["cell_id"] = geo_df.index

    table = TableModel.parse(
        table,
        region_key="region",
        region="polygons",
        instance_key="cell_id",
    )

    sdata.add_shapes("polygons", geo_df, overwrite=True)

    if sdata.table is not None:
        del sdata.table

    sdata.table = table
