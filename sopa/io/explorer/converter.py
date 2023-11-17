import json
import logging
from pathlib import Path

import geopandas as gpd
from spatialdata import SpatialData

from ..._sdata import get_boundaries, get_element, get_spatial_image, to_intrinsic
from . import (
    write_cell_categories,
    write_gene_counts,
    write_image,
    write_polygons,
    write_transcripts,
)
from ._constants import FileNames, experiment_dict
from .utils import explorer_file_path

log = logging.getLogger(__name__)


def _check_explorer_directory(path: Path):
    assert (
        not path.exists() or path.is_dir()
    ), f"A path to an existing file was provided. It should be a path to a directory."
    path.mkdir(parents=True, exist_ok=True)


def write_explorer(
    path: str,
    sdata: SpatialData,
    image_key: str | None = None,
    shapes_key: str | None = None,
    points_key: str | None = None,
    gene_column: str | None = None,
    layer: str | None = None,
    polygon_max_vertices: int = 13,
    lazy: bool = True,
    ram_threshold_gb: int | None = 4,
    save_image_mode: int = 1,
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
        lazy: If `True`, will not load the full images in memory (except if the image memory is below `ram_threshold_gb`).
        ram_threshold_gb: Threshold (in gygabytes) from which image can be loaded in memory. If `None`, the image is never loaded in memory.
        save_image_mode: `1` is normal mode. `0` doesn't save the image. `2` saves **only** the image.
    """
    path: Path = Path(path)
    _check_explorer_directory(path)

    image_key, image = get_spatial_image(sdata, image_key, return_key=True)

    if save_image_mode == 2:
        log.info(f"{save_image_mode=} (only the image will be saved)")
        write_image(path, image, lazy=lazy, ram_threshold_gb=ram_threshold_gb)
        return

    ### Saving cell categories and gene counts
    if sdata.table is not None:
        adata = sdata.table

        shapes_key = adata.uns["spatialdata_attrs"]["region"]
        geo_df = sdata[shapes_key]

        write_gene_counts(path, adata, layer=layer)
        write_cell_categories(path, adata)

    ### Saving cell boundaries
    if shapes_key is None:
        shapes_key, geo_df = get_boundaries(sdata, return_key=True, warn=True)
    else:
        geo_df = sdata[shapes_key]

    if geo_df is not None:
        geo_df = to_intrinsic(sdata, geo_df, image_key)
        write_polygons(path, geo_df.geometry, polygon_max_vertices)

    ### Saving transcripts
    df = get_element(sdata, "points", points_key)
    if df is not None:
        if gene_column is not None:
            df = to_intrinsic(sdata, df, image_key)
            write_transcripts(path, df, gene_column)
        else:
            log.warn("The argument 'gene_column' has to be provided to save the transcripts")

    ### Saving image
    if save_image_mode:
        write_image(path, image, lazy=lazy, ram_threshold_gb=ram_threshold_gb)
    else:
        log.info(f"{save_image_mode=} (the image will not be saved)")

    ### Saving experiment.xenium file
    write_metadata(path, image_key, shapes_key, _get_n_obs(sdata, geo_df))

    log.info(f"Saved all files in the following directory: {path}")
    log.info(f"You can open the experiment with 'open {path / FileNames.METADATA}'")


def _get_n_obs(sdata: SpatialData, geo_df: gpd.GeoDataFrame) -> int:
    if sdata.table is not None:
        return sdata.table.n_obs
    return len(geo_df) if geo_df is not None else 0


def write_metadata(
    path: str, image_key: str = "NA", shapes_key: str = "NA", n_obs: int = 0, is_dir: bool = True
):
    path = explorer_file_path(path, FileNames.METADATA, is_dir)

    with open(path, "w") as f:
        metadata = experiment_dict(image_key, shapes_key, n_obs)
        json.dump(metadata, f, indent=4)
