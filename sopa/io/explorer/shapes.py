from __future__ import annotations

import logging
from math import ceil
from pathlib import Path
from typing import Iterable

import numpy as np
import zarr
from shapely.geometry import Polygon

from ._constants import FileNames, cell_summary_attrs, group_attrs
from .utils import explorer_file_path

log = logging.getLogger(__name__)

TOLERANCE_STEP = 0.5


def pad_polygon(
    polygon: Polygon, max_vertices: int, tolerance: float = TOLERANCE_STEP
) -> np.ndarray:
    """Transform the polygon to have the desired number of vertices

    Args:
        polygon: A `shapely` polygon
        max_vertices: The desired number of vertices
        tolerance: The step of tolerance used for simplification. At each step, we increase the tolerance of this value until the polygon is simplified enough.

    Returns:
        A 2D array representing the polygon vertices
    """
    n_vertices = len(polygon.exterior.coords)
    assert n_vertices >= 3

    coords = polygon.exterior.coords._coords

    if n_vertices == max_vertices:
        return coords.flatten()

    if n_vertices < max_vertices:
        return np.pad(coords, ((0, max_vertices - n_vertices), (0, 0)), mode="edge").flatten()

    # TODO: improve it: how to choose the right tolerance?
    polygon = polygon.simplify(tolerance=tolerance)
    return pad_polygon(polygon, max_vertices, tolerance + TOLERANCE_STEP)


def write_polygons(
    path: Path,
    polygons: Iterable[Polygon],
    max_vertices: int,
    is_dir: bool = True,
    pixel_size: float = 0.2125,
) -> None:
    """Write a `cells.zarr.zip` file containing the cell polygonal boundaries

    Args:
        path: Path to the Xenium Explorer directory where the transcript file will be written
        polygons: A list of `shapely` polygons to be written
        max_vertices: The number of vertices per polygon (they will be transformed to have the right number of vertices)
        is_dir: If `False`, then `path` is a path to a single file, not to the Xenium Explorer directory.
        pixel_size: Number of microns in a pixel. Invalid value can lead to inconsistent scales in the Explorer.
    """
    path = explorer_file_path(path, FileNames.SHAPES, is_dir)

    log.info(f"Writing {len(polygons)} cell polygons")
    coordinates = np.stack([pad_polygon(p, max_vertices) for p in polygons])
    coordinates *= pixel_size

    num_cells = len(coordinates)
    cells_fourth = ceil(num_cells / 4)
    cells_half = ceil(num_cells / 2)

    GROUP_ATTRS = group_attrs()
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
        cell_id[:, 0] = np.arange(num_cells)
        g.array("cell_id", cell_id, dtype="uint32", chunks=(cells_half, 1))

        cell_summary = np.zeros((num_cells, 7))
        cell_summary[:, 2] = [p.area for p in polygons]
        g.array(
            "cell_summary",
            cell_summary,
            dtype="float64",
            chunks=(num_cells, 1),
        )
        g["cell_summary"].attrs.put(cell_summary_attrs())

        g.array(
            "polygon_num_vertices",
            np.full((2, num_cells), n_vertices),
            dtype="int32",
            chunks=(1, cells_half),
        )

        g.array(
            "seg_mask_value",
            np.arange(num_cells),
            dtype="uint32",
            chunks=(cells_half,),
        )
