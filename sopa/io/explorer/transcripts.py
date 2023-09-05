import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

from ...utils.tiling import Tiles2D


def subsample_indices(n_samples, factor: int = 4):
    n_sub = n_samples // factor
    return np.random.choice(n_samples, n_sub, replace=False)


MAX_LEVELS = 15
GRID_SIZE = 250
QUALITY_SCORE = 40


def write_transcripts(path: Path, df: pd.DataFrame):
    num_transcripts = len(df)
    df["gene"] = df["gene"].astype("category")

    location = df[["x", "y"]]
    location = np.concatenate([location, np.zeros((num_transcripts, 1))], axis=1)

    xmin, xmax = location[:, 0].min(), location[:, 0].max()
    ymin, ymax = location[:, 1].min(), location[:, 1].max()

    gene_names = list(df["gene"].cat.categories)
    num_genes = len(gene_names)

    codeword_gene_mapping = list(range(num_genes))

    valid = np.ones((num_transcripts, 1))
    uuid = np.stack(
        [np.arange(num_transcripts), np.full(num_transcripts, 65535)], axis=1
    )
    gene_identity = df["gene"].cat.codes.values[:, None]
    codeword_identity = np.stack(
        [gene_identity[:, 0], np.full(num_transcripts, 65535)], axis=1
    )
    status = np.zeros((num_transcripts, 1))
    quality_score = np.full((num_transcripts, 1), QUALITY_SCORE)

    ATTRS = {
        "codeword_count": num_genes,
        "codeword_gene_mapping": codeword_gene_mapping,
        "codeword_gene_names": gene_names,
        "gene_names": gene_names,
        "gene_index_map": {
            name: index for name, index in zip(gene_names, codeword_gene_mapping)
        },
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
        "grid_size": [GRID_SIZE],
        "grid_array_shapes": [],
        "grid_number_objects": [],
        "grid_keys": [],
    }

    if path.exists():
        path.unlink()

    with zarr.ZipStore(path, mode="w") as store:
        g = zarr.group(store=store)
        g.attrs.put(ATTRS)

        grids = g.create_group("grids")

        for level in range(MAX_LEVELS):
            level_group = grids.create_group(level)

            tile_size = GRID_SIZE * 2**level
            tiles = Tiles2D(xmin, xmax, ymin, ymax, tile_size, 0)

            print(f"Level {level}: {len(tiles)} n_tiles, {len(location)} transcripts")

            indices = tiles.coords_to_indices(location[:, :2])
            tiles_str_indices = np.apply_along_axis(
                lambda l: ",".join(l), 1, indices.astype(str)
            )

            GRIDS_ATTRS["grid_array_shapes"].append([{}] * len(tiles))
            GRIDS_ATTRS["grid_number_objects"].append(len(location))
            GRIDS_ATTRS["grid_keys"].append([])

            for tx in range(tiles.tile_x.count):
                for ty in range(tiles.tile_y.count):
                    str_index = f"{tx},{ty}"
                    loc = np.where(tiles_str_indices == str_index)[0]

                    n_points_tile = len(loc)
                    chunks = (n_points_tile, 1)

                    if n_points_tile == 0:
                        continue

                    GRIDS_ATTRS["grid_keys"][-1].append(str_index)

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

            if len(tiles) == 1:
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

        grids.attrs.put(GRIDS_ATTRS)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        required=True,
        help="Path to the zarr.zip file to be created",
    )
    parser.add_argument(
        "-d",
        "--data",
        type=str,
        required=True,
        help="Path to the pandas transcript file",
    )

    args = parser.parse_args()
    write_transcripts(args.path, pd.read_csv(args.data))
