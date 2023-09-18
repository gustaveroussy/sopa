import argparse

import cv2
import geopandas as gpd
import numpy as np
import spatialdata
import xarray as xr
from cellpose import models
from shapely import affinity
from shapely.geometry import Polygon
from spatialdata.models import ShapesModel
from tqdm import tqdm

from ..utils.tiling import Tiles2D
from .utils import smooth


def extract_polygons(mask: np.ndarray) -> list[Polygon]:
    polys = []

    for cell_id in range(1, mask.max() + 1):
        mask_id = (mask == cell_id).astype("uint8")
        contours, _ = cv2.findContours(mask_id, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)

        polys_ = [smooth(Polygon(c[:, 0, :])) for c in contours if c.shape[0] >= 4]
        polys_ = [p for p in polys_ if not p.is_empty]

        assert len(polys_) <= 1
        polys.extend(polys_)

    return polys


def patch_coordinates(
    xarr: xr.DataArray,
    model: models.CellposeModel,
    x_bounds: list[int],
    y_bounds: list[int],
    c: str,
) -> np.ndarray:
    patch = xarr.sel(
        c=c,
        x=slice(*x_bounds),
        y=slice(*y_bounds),
    ).values

    masks, *_ = model.eval(
        [patch],
        diameter=15,
        channels=[[0, 0]],
        flow_threshold=2,
        cellprob_threshold=-6,
    )
    polygons = extract_polygons(masks[0])

    print(f"Extracted {len(polygons)} polygons")
    return [affinity.translate(p, x_bounds[0], y_bounds[0]) for p in polygons]


def main(args):
    sdata = spatialdata.read_zarr(args.path)
    model = models.Cellpose(model_type="cyto2")

    image = next(iter(sdata.images.values()))

    tiles = Tiles2D(0, len(image.coords["x"]), 0, len(image.coords["y"]), args.width)

    polygons_list = [
        patch_coordinates(image, model, x_bounds, y_bounds, args.dapi)
        for x_bounds, y_bounds in tqdm(tiles[:2])
    ]

    geo_df = gpd.GeoDataFrame({"geometry": [p for l in polygons_list for p in l]})
    geo_df = ShapesModel.parse(geo_df, transformations=image.transform)
    sdata.add_shapes("polygons", geo_df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        required=True,
        help="Path to the zarr spatialdata object",
    )
    parser.add_argument(
        "-d",
        "--dapi",
        type=str,
        required=True,
        help="DAPI channel name",
    )
    parser.add_argument(
        "-w",
        "--width",
        type=int,
        default=5000,
        help="Tile width",
    )

    main(parser.parse_args())
