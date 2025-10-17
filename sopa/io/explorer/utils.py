import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from ..._constants import SopaAttrs
from ...utils import add_spatial_element, get_spatial_element

log = logging.getLogger(__name__)


def explorer_file_path(path: str | Path, filename: str, is_dir: bool):
    path: Path = Path(path)

    if is_dir:
        path = path / filename

    return path


def is_valid_explorer_id(cell_id: str) -> bool:
    """Check if a cell ID is a valid Xenium Explorer ID"""
    if not hasattr(cell_id, "__len__"):
        return False
    if len(cell_id) == 10:
        return cell_id[:-2].isalpha() and cell_id[-2] == "-"
    if len(cell_id) == 8:
        return cell_id.isalpha()
    return False


def int_cell_id(explorer_cell_id: str | pd.Index) -> int | pd.Index:
    """Transforms an alphabetical cell id from the Xenium Explorer to an integer ID

    E.g., int_cell_id('aaaachba-1') = 10000

    Args:
        explorer_cell_id: An alphabetical cell ID or a pandas Index of many explorer cell IDs

    Returns:
        An integer or a pandas Index of integers representing cell IDs as indices"""
    if isinstance(explorer_cell_id, pd.Index):
        return explorer_cell_id.map(int_cell_id)

    assert isinstance(explorer_cell_id, str), "The cell ID must be a string or a pandas Index of strings"
    assert is_valid_explorer_id(explorer_cell_id), "The cell ID must be a valid Xenium Explorer ID"

    code = explorer_cell_id[:-2] if explorer_cell_id[-2] == "-" else explorer_cell_id
    coefs = [ord(c) - 97 for c in code][::-1]
    return sum(value * 16**i for i, value in enumerate(coefs))


def str_cell_id(cell_id: int | pd.Index) -> str | pd.Index:
    """Transforms an integer cell ID into an Xenium Explorer alphabetical cell id

    E.g., str_cell_id(10000) = 'aaaachba-1'

    Args:
        cell_id: An integer or a pandas Index of integers representing cell indices

    Returns:
        A string or a pandas Index of strings representing cell IDs in the Xenium Explorer format
    """
    if isinstance(cell_id, pd.Index):
        return cell_id.map(str_cell_id)

    assert isinstance(cell_id, int), "The cell ID must be an integer or a pandas Index of integers"

    coefs = []
    for _ in range(8):
        cell_id, coef = divmod(cell_id, 16)
        coefs.append(coef)
    return "".join([chr(97 + coef) for coef in coefs][::-1]) + "-1"


def _selection_to_polygon(df, pixel_size):
    return Polygon(df[["X", "Y"]].values / pixel_size)


def xenium_explorer_selection(path: str | Path, pixel_size: float = 0.2125, return_list: bool = False) -> Polygon:
    df = pd.read_csv(path, skiprows=2)

    if "Selection" not in df:
        polygon = _selection_to_polygon(df, pixel_size)
        return [polygon] if return_list else polygon

    return [_selection_to_polygon(sub_df, pixel_size) for _, sub_df in df.groupby("Selection")]


def add_explorer_selection(
    sdata: SpatialData,
    path: str,
    key_added: str = "explorer_selection",
    image_key: str | None = None,
    pixel_size: float = 0.2125,
):
    """After saving a selection on the Xenium Explorer, it will add all polygons inside `sdata.shapes[shapes_key]`

    Args:
        sdata: A `SpatialData` object
        path: The path to the `coordinates.csv` selection file
        key_added: The name to provide to the selection as shapes
        image_key: The original image name
        pixel_size: Number of microns in a pixel. It must be the same value as the one used in `sopa.io.write`
    """
    polys = xenium_explorer_selection(path, pixel_size=pixel_size, return_list=True)
    image = get_spatial_element(sdata.images, key=image_key or sdata.attrs.get(SopaAttrs.CELL_SEGMENTATION))

    transformations = get_transformation(image, get_all=True).copy()

    geo_df = ShapesModel.parse(gpd.GeoDataFrame(geometry=polys), transformations=transformations)
    add_spatial_element(sdata, key_added, geo_df)
