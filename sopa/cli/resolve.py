import typer

app_resolve = typer.Typer()
option = typer.Option()


@app_resolve.command()
def cellpose(
    sdata_path: str,
    patch_dir: str = option,
    expand_radius: float = typer.Option(default=0),
):
    import spatialdata

    from sopa._sdata import get_key
    from sopa.segmentation import StainingSegmentation, shapes
    from sopa.segmentation.cellpose.update import add_shapes

    sdata = spatialdata.read_zarr(sdata_path)

    image_key = get_key(sdata, "images")

    cells = StainingSegmentation.read_patches_cells(patch_dir)
    cells = shapes.solve_conflicts(cells)
    cells = shapes.expand(cells, expand_radius)

    add_shapes(sdata, cells, image_key)


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
    from spatialdata.transformations import get_transformation

    from sopa._constants import SopaKeys
    from sopa._sdata import get_intrinsic_cs, get_item, get_key
    from sopa.segmentation.baysor.aggregate import read_all_baysor_patches, resolve

    sdata = spatialdata.read_zarr(sdata_path)

    patches_cells, adatas = read_all_baysor_patches(baysor_dir, min_area, n)
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
