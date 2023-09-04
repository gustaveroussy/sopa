from math import ceil
from pathlib import Path

import numpy as np
import zarr

CELLS_SUMMARY_ATTRS = {
    "column_descriptions": [
        "Cell centroid in X",
        "Cell centroid in Y",
        "Cell area",
        "Nucleus centroid in X",
        "Nucleus centroid in Y",
        "Nucleus area",
        "z_level",
    ],
    "column_names": [
        "cell_centroid_x",
        "cell_centroid_y",
        "cell_area",
        "nucleus_centroid_x",
        "nucleus_centroid_y",
        "nucleus_area",
        "z_level",
    ],
}

GROUP_ATTRS = {
    "major_version": 5,
    "minor_version": 0,
    "name": "CellSegmentationDataset",
    "polygon_set_descriptions": [
        "DAPI-based nuclei segmentation",
        "Expansion of nuclei boundaries by 15 \u03bcm",
    ],
    "polygon_set_display_names": ["Nucleus boundaries", "Cell boundaries"],
    "polygon_set_names": ["nucleus", "cell"],
    "spatial_units": "microns",
}


def save_polygons(path: Path, coordinates: np.ndarray) -> None:
    num_cells = len(coordinates)
    cells_fourth = ceil(num_cells / 4)
    cells_half = ceil(num_cells / 2)

    GROUP_ATTRS["number_cells"] = num_cells

    polygon_vertices = np.stack([coordinates, coordinates])
    num_points = polygon_vertices.shape[2]
    n_vertices = num_points // 2

    with zarr.ZipStore(path, mode="w") as store:
        g = zarr.group(store=store)
        g.attrs.put(GROUP_ATTRS)

        g.array(
            "polygon_vertices",
            polygon_vertices,
            dtype="float32",
            chunks=(1, cells_fourth, ceil(num_points / 4)),
        )

        cell_id = np.ones((num_cells, 2))
        cell_id[:, 0] = np.arange(1, num_cells + 1)
        g.array("cell_id", cell_id, dtype="uint32", chunks=(cells_half, 1))

        g.array(
            "cell_summary",
            np.zeros((num_cells, 7)),
            dtype="float64",
            chunks=(num_cells, 1),
        )
        g["cell_summary"].attrs.put(CELLS_SUMMARY_ATTRS)

        g.array(
            "polygon_num_vertices",
            np.full((2, num_cells), n_vertices),
            dtype="int32",
            chunks=(1, cells_half),
        )

        g.array(
            "seg_mask_value",
            np.arange(1, num_cells + 1),
            dtype="uint32",
            chunks=(cells_half,),
        )
