import logging
from math import ceil
from pathlib import Path

import dask.dataframe as dd
import numpy as np
import zarr
from zarr.storage import ZipStore

from ._constants import ExplorerConstants, FileNames
from .utils import explorer_file_path

log = logging.getLogger(__name__)


def subsample_indices(indices: np.ndarray, factor: int = 4) -> np.ndarray:
    return np.random.choice(indices, len(indices) // factor, replace=False)


def write_transcripts(
    path: Path,
    df: dd.DataFrame,
    gene: str = "gene",
    max_levels: int = 15,
    is_dir: bool = True,
    pixel_size: float = 0.2125,
):
    """Write a `transcripts.zarr.zip` file containing pyramidal transcript locations

    Args:
        path: Path to the Xenium Explorer directory where the transcript file will be written
        df: DataFrame representing the transcripts, with `"x"`, `"y"` column required, as well as the `gene` column (see the corresponding argument)
        gene: Column of `df` containing the genes names.
        max_levels: Maximum number of levels in the pyramid.
        is_dir: If `False`, then `path` is a path to a single file, not to the Xenium Explorer directory.
        pixel_size: Number of microns in a pixel. Invalid value can lead to inconsistent scales in the Explorer.
    """
    path = explorer_file_path(path, FileNames.POINTS, is_dir)

    # TODO: make everything using dask instead of pandas
    df = df.compute()

    num_transcripts = len(df)
    grid_size = ExplorerConstants.GRID_SIZE / ExplorerConstants.PIXELS_TO_MICRONS * pixel_size
    df[gene] = df[gene].astype("category")

    xy: np.ndarray = (df[["x", "y"]] * pixel_size).values

    if xy.min() < 0:
        log.warning("Some transcripts are located outside of the image (pixels < 0)")
    log.info(f"Writing {len(df)} transcripts")

    xmax, ymax = xy.max(axis=0)

    gene_names = list(df[gene].cat.categories)
    num_genes = len(gene_names)

    codeword_gene_mapping = list(range(num_genes))

    valid = np.ones((num_transcripts, 1), dtype=np.uint8)
    arange = np.arange(num_transcripts, dtype=np.uint32)
    uuid = np.stack([arange, np.full(num_transcripts, 65535, dtype=np.uint32)], axis=1)
    transcript_id = np.stack([arange, np.full(num_transcripts, 65535, dtype=np.uint32)], axis=1)
    gene_identity = df[gene].cat.codes.values[:, None].astype(np.uint16)
    codeword_identity = np.stack([gene_identity[:, 0], np.full(num_transcripts, 65535, dtype=np.uint16)], axis=1)
    status = np.zeros((num_transcripts, 1), dtype=np.uint8)
    quality_score = np.full((num_transcripts, 1), ExplorerConstants.QUALITY_SCORE, dtype=np.float32)

    ATTRS = {
        "codeword_count": num_genes,
        "codeword_gene_mapping": codeword_gene_mapping,
        "codeword_gene_names": gene_names,
        "gene_names": gene_names,
        "gene_index_map": dict(zip(gene_names, codeword_gene_mapping)),
        "number_genes": num_genes,
        "spatial_units": "micron",
        "coordinate_space": "refined-final_global_micron",
        "major_version": 4,
        "minor_version": 1,
        "name": "RnaDataset",
        "number_rnas": num_transcripts,
        "dataset_uuid": "unique-id-test",
        "data_format": 0,
    }

    GRIDS_ATTRS = {
        "grid_key_names": ["grid_x_loc", "grid_y_loc"],
        "grid_zip": False,
        "grid_size": [grid_size],
        "grid_array_shapes": [],
        "grid_number_objects": [],
        "grid_keys": [],
    }

    subsampling_locs = {0: np.arange(num_transcripts)}

    for level in range(max_levels):
        tile_size = grid_size * 2**level
        level_xy = xy[subsampling_locs[level]]

        indices = np.floor(level_xy / tile_size).clip(0).astype(int)
        tiles_str_indices = np.array([f"{tx},{ty}" for (tx, ty) in indices])

        GRIDS_ATTRS["grid_array_shapes"].append([])
        GRIDS_ATTRS["grid_number_objects"].append([])
        GRIDS_ATTRS["grid_keys"].append([])

        n_tiles_x, n_tiles_y = ceil(xmax / tile_size), ceil(ymax / tile_size)

        for tx in range(n_tiles_x):
            for ty in range(n_tiles_y):
                str_index = f"{tx},{ty}"
                loc = np.where(tiles_str_indices == str_index)[0]

                n_points_tile = len(loc)

                if n_points_tile == 0:
                    continue

                GRIDS_ATTRS["grid_array_shapes"][-1].append({})
                GRIDS_ATTRS["grid_keys"][-1].append(str_index)
                GRIDS_ATTRS["grid_number_objects"][-1].append(n_points_tile)

        if n_tiles_x * n_tiles_y == 1:
            GRIDS_ATTRS["number_levels"] = level + 1
            break

        if level + 1 < max_levels:
            subsampling_locs[level + 1] = subsample_indices(subsampling_locs[level])

    with ZipStore(path, mode="w") as store:
        g = zarr.group(store=store, zarr_format=2, attributes=ATTRS)

        grids = g.create_group("grids", attributes=GRIDS_ATTRS)

        for level, level_locs in subsampling_locs.items():
            log.info(f"   > Level {level}: {len(level_locs)} transcripts")
            level_group = grids.create_group(str(level))

            tile_size = grid_size * 2**level

            indices = np.floor(xy[level_locs] / tile_size).clip(0).astype(int)
            tiles_str_indices = np.array([f"{tx},{ty}" for (tx, ty) in indices])

            n_tiles_x, n_tiles_y = ceil(xmax / tile_size), ceil(ymax / tile_size)

            for tx in range(n_tiles_x):
                for ty in range(n_tiles_y):
                    str_index = f"{tx},{ty}"
                    loc = np.where(tiles_str_indices == str_index)[0]

                    n_points_tile = len(loc)
                    chunks = (n_points_tile, 1)

                    if n_points_tile == 0:
                        continue

                    location = np.concatenate([xy[level_locs][loc], np.zeros((len(loc), 1))], axis=1).astype(np.float32)

                    tile_group = level_group.create_group(str_index)
                    tile_group.create_array("valid", data=valid[level_locs][loc], chunks=chunks)
                    tile_group.create_array("status", data=status[level_locs][loc], chunks=chunks)
                    tile_group.create_array("location", data=location, chunks=chunks)
                    tile_group.create_array("gene_identity", data=gene_identity[level_locs][loc], chunks=chunks)
                    tile_group.create_array("quality_score", data=quality_score[level_locs][loc], chunks=chunks)
                    tile_group.create_array("codeword_identity", data=codeword_identity[level_locs][loc], chunks=chunks)
                    tile_group.create_array("uuid", data=uuid[level_locs][loc], chunks=chunks)
                    tile_group.create_array("id", data=transcript_id[level_locs][loc], chunks=chunks)
