from __future__ import annotations

import logging
from math import ceil
from pathlib import Path

import dask.dataframe as dd
import numpy as np
import zarr

from ._constants import ExplorerConstants, FileNames
from .utils import explorer_file_path

log = logging.getLogger(__name__)


def subsample_indices(n_samples, factor: int = 4):
    n_sub = n_samples // factor
    return np.random.choice(n_samples, n_sub, replace=False)


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

    location = df[["x", "y"]]
    location *= pixel_size
    location = np.concatenate([location, np.zeros((num_transcripts, 1))], axis=1)

    if location.min() < 0:
        log.warn("Some transcripts are located outside of the image (pixels < 0)")
    log.info(f"Writing {len(df)} transcripts")

    xmax, ymax = location[:, :2].max(axis=0)

    gene_names = list(df[gene].cat.categories)
    num_genes = len(gene_names)

    codeword_gene_mapping = list(range(num_genes))

    valid = np.ones((num_transcripts, 1))
    uuid = np.stack([np.arange(num_transcripts), np.full(num_transcripts, 65535)], axis=1)
    transcript_id = np.stack([np.arange(num_transcripts), np.full(num_transcripts, 65535)], axis=1)
    gene_identity = df[gene].cat.codes.values[:, None]
    codeword_identity = np.stack([gene_identity[:, 0], np.full(num_transcripts, 65535)], axis=1)
    status = np.zeros((num_transcripts, 1))
    quality_score = np.full((num_transcripts, 1), ExplorerConstants.QUALITY_SCORE)

    ATTRS = {
        "codeword_count": num_genes,
        "codeword_gene_mapping": codeword_gene_mapping,
        "codeword_gene_names": gene_names,
        "gene_names": gene_names,
        "gene_index_map": {name: index for name, index in zip(gene_names, codeword_gene_mapping)},
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

    with zarr.ZipStore(path, mode="w") as store:
        g = zarr.group(store=store)
        g.attrs.put(ATTRS)

        grids = g.create_group("grids")

        for level in range(max_levels):
            log.info(f"   > Level {level}: {len(location)} transcripts")
            level_group = grids.create_group(level)

            tile_size = grid_size * 2**level

            indices = np.floor(location[:, :2] / tile_size).clip(0).astype(int)
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
                    chunks = (n_points_tile, 1)

                    if n_points_tile == 0:
                        continue

                    GRIDS_ATTRS["grid_array_shapes"][-1].append({})
                    GRIDS_ATTRS["grid_keys"][-1].append(str_index)
                    GRIDS_ATTRS["grid_number_objects"][-1].append(n_points_tile)

                    tile_group = level_group.create_group(str_index)
                    tile_group.array(
                        "valid",
                        valid[loc],
                        dtype="uint8",
                        chunks=chunks,
                    )
                    tile_group.array(
                        "status",
                        status[loc],
                        dtype="uint8",
                        chunks=chunks,
                    )
                    tile_group.array(
                        "location",
                        location[loc],
                        dtype="float32",
                        chunks=chunks,
                    )
                    tile_group.array(
                        "gene_identity",
                        gene_identity[loc],
                        dtype="uint16",
                        chunks=chunks,
                    )
                    tile_group.array(
                        "quality_score",
                        quality_score[loc],
                        dtype="float32",
                        chunks=chunks,
                    )
                    tile_group.array(
                        "codeword_identity",
                        codeword_identity[loc],
                        dtype="uint16",
                        chunks=chunks,
                    )
                    tile_group.array(
                        "uuid",
                        uuid[loc],
                        dtype="uint32",
                        chunks=chunks,
                    )
                    tile_group.array(
                        "id",
                        transcript_id[loc],
                        dtype="uint32",
                        chunks=chunks,
                    )

            if n_tiles_x * n_tiles_y == 1:
                GRIDS_ATTRS["number_levels"] = level + 1
                break

            sub_indices = subsample_indices(len(location))

            location = location[sub_indices]
            valid = valid[sub_indices]
            status = status[sub_indices]
            gene_identity = gene_identity[sub_indices]
            quality_score = quality_score[sub_indices]
            codeword_identity = codeword_identity[sub_indices]
            uuid = uuid[sub_indices]
            transcript_id = transcript_id[sub_indices]

        grids.attrs.put(GRIDS_ATTRS)
