import geopandas as gpd
import numpy as np
import pytest
from shapely.affinity import translate
from shapely.geometry import Polygon

from sopa.segmentation import solve_conflicts
from sopa.segmentation.shapes import vectorize


@pytest.fixture
def mask() -> np.ndarray:
    return np.load("tests/mask_example.npy")


@pytest.fixture
def cells(mask: np.ndarray) -> gpd.GeoDataFrame:
    return vectorize(mask)


def test_keep_all_cells(mask: np.ndarray, cells: gpd.GeoDataFrame):
    assert len(cells) == mask.max()


def test_all_polygons(cells: gpd.GeoDataFrame):
    assert all(isinstance(cell, Polygon) for cell in cells.geometry)


def test_solve_conflict(cells: gpd.GeoDataFrame):
    other_cells = [translate(cell, 2, 3) for cell in cells.geometry]

    res = solve_conflicts(list(cells.geometry) + other_cells)
    assert all(isinstance(cell, Polygon) for cell in res.geometry)


def test_solve_conflict_patches(cells: gpd.GeoDataFrame):
    tile_overlap = 50
    other_cells = [translate(cell, cells.total_bounds[2] - tile_overlap, 0) for cell in cells.geometry]

    res = solve_conflicts(
        cells=list(cells.geometry) + other_cells,
        patch_indices=np.array([0] * len(cells) + [1] * len(other_cells)),
        patch_centroids={
            0: (cells.total_bounds[2] / 2, cells.total_bounds[3] / 2),
            1: (cells.total_bounds[2] / 2 + (cells.total_bounds[2] - tile_overlap), cells.total_bounds[3] / 2),
        },
    )
    assert all(isinstance(cell, Polygon) for cell in res.geometry)
    assert len(res) < len(cells) + len(other_cells)
