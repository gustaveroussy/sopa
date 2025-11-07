import geopandas as gpd
import numpy as np
import pytest
from shapely.affinity import translate
from shapely.geometry import Polygon

import sopa
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

    res: Polygon = _vectorize_mask(mask, allow_holes=False).geometry[0]
    assert len(res.interiors) == 0

    res: Polygon = _vectorize_mask(mask, allow_holes=True).geometry[0]
    assert len(res.interiors) == 2


def test_simplification_settings():
    mask1 = np.zeros((10, 10), dtype=np.uint8)
    mask1[2:6, 2:6] = 1
    mask1[4:9, 4:9] = 1
    mask1[3, 6] = 1
    mask1[1, 2] = 1

    radius = 30
    center = (radius + 1, radius + 1)
    size = 2 * radius + 3
    mask2 = np.zeros((size, size), dtype=np.uint8)
    y, x = np.ogrid[-center[0] : size - center[0], -center[1] : size - center[1]]
    mask2[x**2 + y**2 <= radius**2] = 1

    for mask in [mask1, mask2]:
        sopa.settings.simplification_tolerance = None

        cell_normal = _vectorize_mask(mask).geometry[0]
        cell_low_simp = vectorize(mask, 0.05).geometry[0]

        sopa.settings.simplification_tolerance = 10
        cell_simple = vectorize(mask).geometry[0]

        sopa.settings.simplification_tolerance = 0
        cell_no_simp = vectorize(mask).geometry[0]

        iou = cell_normal.intersection(cell_no_simp).area / cell_normal.union(cell_no_simp).area
        assert iou > 0.98

        iou = cell_normal.intersection(cell_low_simp).area / cell_normal.union(cell_low_simp).area
        assert iou > 0.98

        iou = cell_normal.intersection(cell_simple).area / cell_normal.union(cell_simple).area
        assert iou < 0.8

    sopa.settings.simplification_tolerance = None
