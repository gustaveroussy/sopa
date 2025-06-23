import geopandas as gpd
import numpy as np
import pytest
from shapely.affinity import translate
from shapely.geometry import Polygon

from sopa.segmentation import solve_conflicts
from sopa.shapes import _vectorize_mask, vectorize


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


def test_vectorize_with_holes():
    mask = np.zeros((20, 20), dtype=np.uint8)
    mask[2:18, 2:18] = 1
    mask[3:6, 10:16] = 0  # hole 1
    mask[8:14, 8:14] = 0  # hole 2

    res = _vectorize_mask(mask, allow_holes=False).geometry[0]
    assert len(res.geoms[0].interiors) == 0

    res = _vectorize_mask(mask, allow_holes=True).geometry[0]
    assert len(res.geoms[0].interiors) == 2
