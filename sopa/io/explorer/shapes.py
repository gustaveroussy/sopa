from math import ceil
from pathlib import Path
from typing import Iterable

import numpy as np
import zarr
from shapely.geometry import Polygon

from ._constants import ExplorerConstants, cell_summary_attrs, group_attrs


def pad_polygon(polygon: Polygon, max_vertices: int, tolerance: float = 1) -> np.ndarray:
    n_vertices = len(polygon.exterior.coords)
    assert n_vertices >= 3

    coords = polygon.exterior.coords._coords

    if n_vertices == max_vertices:
        return coords.flatten()

    if n_vertices < max_vertices:
        return np.pad(coords, ((0, max_vertices - n_vertices), (0, 0)), mode="edge").flatten()

    # TODO: improve it: how to choose the right tolerance?
    polygon = polygon.simplify(tolerance=tolerance)
    return pad_polygon(polygon, max_vertices, tolerance + 1)


def write_polygons(path: Path, polygons: Iterable[Polygon], max_vertices: int) -> None:
    print(f"Writing {len(polygons)} cell polygons")
    coordinates = np.stack([pad_polygon(p, max_vertices) for p in polygons])
    coordinates /= ExplorerConstants.MICRONS_TO_PIXELS

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
        cell_id[:, 0] = np.arange(1, num_cells + 1)
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
            np.arange(1, num_cells + 1),
            dtype="uint32",
            chunks=(cells_half,),
        )
