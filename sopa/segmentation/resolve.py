import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import shapely.affinity
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from tqdm import tqdm

from ..utils import to_intrinsic
from .shapes import _ensure_polygon

log = logging.getLogger(__name__)


def solve_conflicts(
    cells: list[Polygon] | gpd.GeoDataFrame,
    threshold: float = 0.5,
    patch_indices: np.ndarray | None = None,
    return_indices: bool = False,
) -> gpd.GeoDataFrame | tuple[gpd.GeoDataFrame, np.ndarray]:
    """Resolve segmentation conflicts (i.e. overlap) after running segmentation on patches

    Args:
        cells: List of cell polygons
        threshold: When two cells are overlapping, we look at the area of intersection over the area of the smallest cell. If this value is higher than the `threshold`, the cells are merged
        patch_indices: Patch from which each cell belongs.
        return_indices: If `True`, returns also the cells indices. Merged cells have an index of -1.

    Returns:
        Array of resolved cells polygons. If `return_indices`, it also returns an array of cell indices.
    """
    cells = list(cells.geometry) if isinstance(cells, gpd.GeoDataFrame) else list(cells)
    n_cells = len(cells)
    resolved_indices = np.arange(n_cells)

    assert n_cells > 0, "No cells was segmented, cannot continue"

    tree = shapely.STRtree(cells)
    conflicts = tree.query(cells, predicate="intersects")

    if patch_indices is not None:
        conflicts = conflicts[:, patch_indices[conflicts[0]] != patch_indices[conflicts[1]]].T
    else:
        conflicts = conflicts[:, conflicts[0] != conflicts[1]].T

    for i1, i2 in tqdm(conflicts, desc="Resolving conflicts"):
        resolved_i1: int = resolved_indices[i1]
        resolved_i2: int = resolved_indices[i2]
        cell1, cell2 = cells[resolved_i1], cells[resolved_i2]

        intersection = cell1.intersection(cell2).area
        if intersection >= threshold * min(cell1.area, cell2.area):
            cell = _ensure_polygon(cell1.union(cell2))
            assert not cell.is_empty, "Merged cell is empty"

            resolved_indices[np.isin(resolved_indices, [resolved_i1, resolved_i2])] = len(cells)
            cells.append(cell)

    unique_indices = np.unique(resolved_indices)
    unique_cells = gpd.GeoDataFrame(geometry=cells).iloc[unique_indices]

    if return_indices:
        return unique_cells, np.where(unique_indices < n_cells, unique_indices, -1)

    return unique_cells


def combine(
    sdata: SpatialData,
    elements: list[str | gpd.GeoDataFrame],
    key_added: str,
    threshold: float = 0.5,
):
    """Combine multiple segmentation boundaries into a single one.

    Example:
        On the example below, we run Cellpose twice, once for nuclei and once for tumor cells. We then combine the two segmentations into a single one.
        ```python
        import sopa

        sdata = sopa.io.toy_dataset(length=1000)
        sopa.make_image_patches(sdata)

        sopa.segmentation.cellpose(sdata, "DAPI", diameter=35, key_added="nuclei")
        sopa.segmentation.cellpose(sdata, ["DAPI", "CK"], diameter=35, key_added="tumor_cells")

        sopa.segmentation.combine(sdata, ["nuclei", "tumor_cells"], key_added="combined_cells")
        ```

    Args:
        sdata: A `SpatialData` object.
        elements: List of name of the keys in `sdata.shapes` to be combined (or directly a list of `GeoDataFrame`).
        key_added: The name of the new key to be added to `sdata.shapes`.
        threshold: When two cells are overlapping, we look at the area of intersection over the area of the smallest cell. If this value is higher than the `threshold`, the cells are merged
    """
    assert len(elements) > 1, "At least two elements must be provided to combine"

    elements: list[gpd.GeoDataFrame] = [
        element if isinstance(element, gpd.GeoDataFrame) else sdata.shapes[element] for element in elements
    ]

    reference = elements[0]
    intrinsic_elements = [reference] + [to_intrinsic(sdata, element, reference) for element in elements[1:]]

    combined_cells = list(pd.concat([element.geometry for element in intrinsic_elements], axis=0))
    combined_cells = solve_conflicts(combined_cells, threshold=threshold)

    combined_geo_df = ShapesModel.parse(combined_cells, transformations=get_transformation(reference, get_all=True))

    sdata.shapes[key_added] = combined_geo_df
