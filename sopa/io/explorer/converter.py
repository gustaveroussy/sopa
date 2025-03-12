import json
import logging
from pathlib import Path

import geopandas as gpd
from anndata import AnnData
from spatialdata import SpatialData

from ..._constants import ATTRS_KEY, SopaAttrs, SopaKeys
from ...utils import (
    get_boundaries,
    get_feature_key,
    get_spatial_element,
    get_spatial_image,
    to_intrinsic,
)
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
    ), "A path to an existing file was provided. It should be a path to a directory."
    path.mkdir(parents=True, exist_ok=True)


def _should_save(mode: str | None, character: str):
    if mode is None:
        return True

    assert len(mode) and mode[0] in [
        "-",
        "+",
    ], "Mode should be a string that starts with '+' or '-'"

    return character in mode if mode[0] == "+" else character not in mode


def write(
    path: str,
    sdata: SpatialData,
    table_key: str = SopaKeys.TABLE,
    image_key: str | None = None,
    shapes_key: str | None = None,
    points_key: str | None = None,
    gene_column: str | None = None,
    pixel_size: float = 0.2125,
    layer: str | None = None,
    polygon_max_vertices: int = 13,
    lazy: bool = True,
    ram_threshold_gb: int | None = 4,
    mode: str = None,
    save_h5ad: bool = False,
    run_name: str | None = None,
) -> None:
    """
    Transform a SpatialData object into inputs for the Xenium Explorer.
    After running this function, double-click on the `experiment.xenium` file to open it.

    !!! note "Software download"
        Make sure you have the latest version of the [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer)

    Note:
        This function will create up to 7 files, depending on the `SpatialData` object and the arguments:

        - `experiment.xenium` contains some experiment metadata. Double-click on this file to open the Xenium Explorer. This file can also be created with [`write_metadata`](./#sopa.io.explorer.write_metadata).

        - `morphology.ome.tif` is the primary image. This file can also be created with [`write_image`](./#sopa.io.explorer.write_image). Add more images with `align`.

        - `analysis.zarr.zip` contains the cells categories (or clusters), i.e. `adata.obs`. This file can also be created with [`write_cell_categories`](./#sopa.io.explorer.write_cell_categories).

        - `cell_feature_matrix.zarr.zip` contains the cell-by-gene counts. This file can also be created with [`write_gene_counts`](./#sopa.io.explorer.write_gene_counts).

        - `cells.zarr.zip` contains the cells polygon boundaries. This file can also be created with [`write_polygons`](./#sopa.io.explorer.write_polygons).

        - `transcripts.zarr.zip` contains transcripts locations. This file can also be created with [`write_transcripts`](./#sopa.io.explorer.write_transcripts).

        - `adata.h5ad` is the `AnnData` object from the `SpatialData`. This is **not** used by the Explorer, but only saved for convenience.

    Args:
        path: Path to the directory where files will be saved.
        sdata: SpatialData object.
        table_key: Name of the table containing the gene counts or intensities (key of `sdata.tables`). By default, uses `sdata["table"]`.
        image_key: Name of the image of interest (key of `sdata.images`). By default, it will be inferred.
        shapes_key: Name of the cell shapes (key of `sdata.shapes`). By default, it will be inferred from the table.
        points_key: Name of the transcripts (key of `sdata.points`). By default, it will be inferred.
        gene_column: Column name of the points dataframe containing the gene names.
        pixel_size: Number of microns in a pixel. Invalid value can lead to inconsistent scales in the Explorer.
        layer: Layer of the AnnData table where the gene counts are saved. If `None`, uses `table.X`.
        polygon_max_vertices: Maximum number of vertices for the cell polygons.
        lazy: If `True`, will not load the full images in memory (except if the image memory is below `ram_threshold_gb`).
        ram_threshold_gb: Threshold (in gygabytes) from which image can be loaded in memory. If `None`, the image is never loaded in memory.
        mode: string that indicated which files should be created. "-ib" means everything except images and boundaries, while "+tocm" means only transcripts/observations/counts/metadata (each letter corresponds to one explorer file). By default, keeps everything.
        save_h5ad: Whether to save the adata as h5ad in the explorer directory (for convenience only, since h5ad is faster to open than the original .zarr table)
        run_name: Name of the run displayed in the Xenium Explorer. If `None`, uses the `image_key`.
    """
    path: Path = Path(path)
    _check_explorer_directory(path)

    image_key, _ = get_spatial_image(sdata, key=image_key, return_key=True)

    ### Saving table / cell categories / gene counts
    if table_key in sdata.tables:
        adata: AnnData = sdata.tables[table_key]

        _shapes_key = adata.uns[ATTRS_KEY]["region"]
        assert (
            shapes_key is None or _shapes_key == shapes_key
        ), f"Got {shapes_key=}, while the table corresponds to the shapes {_shapes_key}"
        shapes_key = _shapes_key[0] if isinstance(_shapes_key, list) else _shapes_key

        geo_df = sdata[shapes_key]

        if _should_save(mode, "c"):
            write_gene_counts(path, adata, layer=layer)
        if _should_save(mode, "o"):
            write_cell_categories(path, adata)
        if save_h5ad:
            adata.write_h5ad(path / FileNames.H5AD)

    ### Saving cell boundaries
    if shapes_key is None:
        shapes_key, geo_df = get_boundaries(sdata, return_key=True, warn=True)
    else:
        geo_df = sdata[shapes_key]

    if _should_save(mode, "b") and geo_df is not None:
        geo_df = to_intrinsic(sdata, geo_df, image_key)

        if table_key in sdata.tables:
            geo_df = geo_df.loc[adata.obs[adata.uns[ATTRS_KEY]["instance_key"]]]

        write_polygons(path, geo_df.geometry, polygon_max_vertices, pixel_size=pixel_size)

    ### Saving transcripts
    df = None
    if len(sdata.points):
        df = get_spatial_element(sdata.points, key=points_key or sdata.attrs.get(SopaAttrs.TRANSCRIPTS))

    if _should_save(mode, "t") and not _use_symlink(path, sdata, "transcripts*") and df is not None:
        gene_column = gene_column or get_feature_key(df)
        if gene_column is not None:
            df = to_intrinsic(sdata, df, image_key)
            write_transcripts(path, df, gene_column, pixel_size=pixel_size)
        else:
            log.warning("The argument 'gene_column' has to be provided to save the transcripts")

    ### Saving image
    if _should_save(mode, "i") and not _use_symlink(path, sdata, "morphology*"):
        write_image(
            path,
            sdata[image_key],
            lazy=lazy,
            ram_threshold_gb=ram_threshold_gb,
            pixel_size=pixel_size,
        )

    ### Saving experiment.xenium file
    if _should_save(mode, "m"):
        write_metadata(path, run_name or image_key, shapes_key, _get_n_obs(sdata, geo_df, table_key), pixel_size)

    log.info(f"Saved files in the following directory: {path}")
    log.info(f"You can open the experiment with 'open {path / FileNames.METADATA}'")


def _use_symlink(path: Path, sdata: SpatialData, pattern: str) -> bool:
    """Try using the Xenium output files when existing to avoid re-generating large files."""
    if SopaAttrs.XENIUM_OUTPUT_PATH not in sdata.attrs:
        return False

    files = list(Path(sdata.attrs[SopaAttrs.XENIUM_OUTPUT_PATH]).glob(pattern))
    for file in files:
        target = path / file.name

        if target.exists():
            if not target.is_symlink():  # avoid removing non-symlink files
                return False
            target.unlink()

        target.symlink_to(file)
        log.info(f"Created symlink {target} -> {file}")

    return len(files) > 0


def _get_n_obs(sdata: SpatialData, geo_df: gpd.GeoDataFrame, table_key: str) -> int:
    if table_key in sdata.tables:
        return sdata.tables[table_key].n_obs
    return len(geo_df) if geo_df is not None else 0


def write_metadata(
    path: str,
    run_name: str = "NA",
    shapes_key: str = "NA",
    n_obs: int = 0,
    is_dir: bool = True,
    pixel_size: float = 0.2125,
):
    """Create an `experiment.xenium` file that can be open by the Xenium Explorer.

    Note:
        This function alone is not enough to actually open an experiment. You will need at least to wrun `write_image`, or create all the outputs with `write`.

    Args:
        path: Path to the Xenium Explorer directory where the metadata file will be written
        run_name: Key of `SpatialData` object containing the primary image used on the explorer.
        shapes_key: Key of `SpatialData` object containing the boundaries shown on the explorer.
        n_obs: Number of cells
        is_dir: If `False`, then `path` is a path to a single file, not to the Xenium Explorer directory.
        pixel_size: Number of microns in a pixel. Invalid value can lead to inconsistent scales in the Explorer.
    """
    path = explorer_file_path(path, FileNames.METADATA, is_dir)

    with open(path, "w") as f:
        metadata = experiment_dict(run_name, shapes_key, n_obs, pixel_size)
        json.dump(metadata, f, indent=4)
