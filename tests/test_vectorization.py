import geopandas as gpd
import numpy as np
import pytest
from shapely.affinity import translate
from shapely.geometry import Polygon

from sopa.segmentation import solve_conflicts
from sopa.segmentation.shapes import geometrize


@pytest.fixture
def mask() -> np.ndarray:
    return np.load("tests/mask_example.npy")


@pytest.fixture
def cells(mask: np.ndarray) -> gpd.GeoDataFrame:
    return geometrize(mask)


def test_keep_all_cells(mask: np.ndarray, cells: gpd.GeoDataFrame):
    assert len(cells) == mask.max()


def test_all_polygons(cells: gpd.GeoDataFrame):
    assert all(isinstance(cell, Polygon) for cell in cells.geometry)


def test_solve_conflict(cells: gpd.GeoDataFrame):
    other_cells = [translate(cell, 2, 3) for cell in cells.geometry]

    res = solve_conflicts(list(cells.geometry) + other_cells)
    assert all(isinstance(cell, Polygon) for cell in res.geometry)
