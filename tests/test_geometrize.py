import numpy as np
import pytest
from shapely.affinity import translate
from shapely.geometry import Polygon

from sopa.segmentation.shapes import geometrize, solve_conflicts


@pytest.fixture
def mask() -> np.ndarray:
    return np.load("tests/mask_example.npy")


@pytest.fixture
def cells(mask: np.ndarray):
    return geometrize(mask)


def test_keep_all_cells(mask: np.ndarray, cells: list[Polygon]):
    assert len(cells) == mask.max()


def test_all_polygons(cells: list[Polygon]):
    assert all(isinstance(cell, Polygon) for cell in cells)


def test_solve_conflict(cells: list[Polygon]):
    other_cells = [translate(cell, 2, 3) for cell in cells]

    res = solve_conflicts(cells + other_cells)
    assert all(isinstance(cell, Polygon) for cell in res)
