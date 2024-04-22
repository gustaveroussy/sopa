from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely import Polygon
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from ..._sdata import get_element


def explorer_file_path(path: str, filename: str, is_dir: bool):
    path: Path = Path(path)

    if is_dir:
        path = path / filename

    return path


def int_cell_id(explorer_cell_id: str) -> int:
    """Transforms an alphabetical cell id from the Xenium Explorer to an integer ID

    E.g., int_cell_id('aaaachba-1') = 10000"""
    code = explorer_cell_id[:-2] if explorer_cell_id[-2] == "-" else explorer_cell_id
    coefs = [ord(c) - 97 for c in code][::-1]
    return sum(value * 16**i for i, value in enumerate(coefs))


def str_cell_id(cell_id: int) -> str:
    """Transforms an integer cell ID into an Xenium Explorer alphabetical cell id

    E.g., str_cell_id(10000) = 'aaaachba-1'"""
    coefs = []
    for _ in range(8):
        cell_id, coef = divmod(cell_id, 16)
        coefs.append(coef)
    return "".join([chr(97 + coef) for coef in coefs][::-1]) + "-1"


def _selection_to_polygon(df, pixel_size):
    return Polygon(df[["X", "Y"]].values / pixel_size)


def xenium_explorer_selection(
    path: str | Path, pixel_size: float = 0.2125, return_list: bool = False
) -> Polygon:
    df = pd.read_csv(path, skiprows=2)

    if "Selection" not in df:
        polygon = _selection_to_polygon(df, pixel_size)
        return [polygon] if return_list else polygon

    return [_selection_to_polygon(sub_df, pixel_size) for _, sub_df in df.groupby("Selection")]


def add_explorer_selection(
    sdata: SpatialData,
    path: str,
    shapes_key: str,
    image_key: str | None = None,
    pixel_size: float = 0.2125,
):
    """After saving a selection on the Xenium Explorer, it will add all polygons inside `sdata.shapes[shapes_key]`

    Args:
        sdata: A `SpatialData` object
        path: The path to the `coordinates.csv` selection file
        shapes_key: The name to provide to the shapes
        image_key: The original image name
        pixel_size: Number of microns in a pixel. It must be the same value as the one used in `sopa.io.write`
    """
    polys = xenium_explorer_selection(path, pixel_size=pixel_size, return_list=True)
    image = get_element(sdata, "images", image_key)

    transformations = get_transformation(image, get_all=True).copy()

    sdata.shapes[shapes_key] = ShapesModel.parse(
        gpd.GeoDataFrame(geometry=polys), transformations=transformations
    )
