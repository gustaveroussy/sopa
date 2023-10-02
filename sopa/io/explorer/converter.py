import json
from pathlib import Path

from spatialdata import SpatialData

from ...utils.utils import _get_element, _get_key, get_intrinsic_cs
from . import (
    write_cell_categories,
    write_gene_counts,
    write_image,
    write_polygons,
    write_transcripts,
)
from ._constants import FileNames, experiment_dict


def _reorder_instances(sdata: SpatialData, shapes_key: str):
    adata = sdata.table

    instance_key = adata.uns["spatialdata_attrs"]["instance_key"]
    region_key = adata.uns["spatialdata_attrs"]["region_key"]

    adata = adata[adata.obs[region_key] == shapes_key].copy()
    adata.obs.set_index(instance_key, inplace=True)
    adata = adata[sdata.shapes[shapes_key].index].copy()
    return adata


def write_explorer(
    path: str,
    sdata: SpatialData,
    image_key: str | None = None,
    shapes_key: str | None = None,
    points_key: str | None = None,
    gene_column: str | None = None,
    layer: str | None = None,
    polygon_max_vertices: int = 13,
) -> None:
    """
    Transform a SpatialData object into inputs for the Xenium Explorer.
    Currently only images of type MultiscaleSpatialImage are supported.

    Args:
        path: Path to the directory where files will be saved.
        sdata: SpatialData object.
        image_key: Name of the image of interest (key of `sdata.images`).
        shapes_key: Name of the cell shapes (key of `sdata.shapes`).
        points_key: Name of the transcripts (key of `sdata.points`).
        gene_column: Column name of the points dataframe containing the gene names.
        layer: Layer of `sdata.table` where the gene counts are saved. If `None`, uses `sdata.table.X`.
        polygon_max_vertices: Maximum number of vertices for the cell polygons.
    """
    path: Path = Path(path)
    assert (
        not path.exists() or path.is_dir()
    ), f"A path to an existing file was provided. It should be a path to a directory."

    path.mkdir(parents=True, exist_ok=True)

    image_key = _get_key(sdata, "images", image_key)
    assert image_key is not None, "An image is required to convert to the Xenium Explorer inputs"

    if sdata.table is not None:
        adata = _reorder_instances(sdata, shapes_key)

        write_gene_counts(path / FileNames.TABLE, adata, layer)
        write_cell_categories(path / FileNames.CELL_CATEGORIES, adata)

    pixels_cs = get_intrinsic_cs(sdata, image_key)

    shapes_key = _get_key(sdata, "shapes", shapes_key)
    gdf = sdata.shapes[shapes_key]
    if gdf is not None:
        gdf = sdata.transform_element_to_coordinate_system(gdf, pixels_cs)
        write_polygons(path / FileNames.SHAPES, gdf.geometry, polygon_max_vertices)

    df = _get_element(sdata, "points", points_key)
    if df is not None:
        assert (
            gene_column is not None
        ), "The argument 'gene_column' has to be provided to save the transcripts"
        df = sdata.transform_element_to_coordinate_system(df, pixels_cs)
        write_transcripts(path / FileNames.POINTS, df, gene_column)

    write_image(path / FileNames.IMAGE, sdata.images[image_key], image_key)

    n_obs = sdata.table.n_obs if sdata.table is not None else (len(gdf) if gdf is not None else 0)

    EXPERIMENT = experiment_dict(image_key, shapes_key, n_obs)
    with open(path / FileNames.METADATA, "w") as f:
        json.dump(EXPERIMENT, f, indent=4)
