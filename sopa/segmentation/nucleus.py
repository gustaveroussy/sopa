import argparse

import geopandas as gpd
import numpy as np
import spatialdata
import xarray as xr
from cellpose import models
from shapely import affinity
from spatialdata.models import ShapesModel
from tqdm import tqdm

from ..utils.tiling import Tiles2D
from .utils import extract_polygons, solve_conflicts


def cellpose_patch(model_type: str = "cyto2"):
    model = models.Cellpose(model_type=model_type)

    def _(patch: np.ndarray):
        masks, *_ = model.eval(
            [patch],
            diameter=15,
            channels=[[0, 0]],
            flow_threshold=2,
            cellprob_threshold=-6,
        )
        return masks[0]

    return _


def run_patch(
    xarr: xr.DataArray,
    x_bounds: list[int],
    y_bounds: list[int],
    c: str,
    method: callable,
) -> np.ndarray:
    patch = xarr.sel(
        c=c,
        x=slice(*x_bounds),
        y=slice(*y_bounds),
    ).values

    polygons = extract_polygons(method(patch))

    return [affinity.translate(p, x_bounds[0], y_bounds[0]) for p in polygons]


def main(args):
    sdata = spatialdata.read_zarr(args.path)

    image = next(iter(sdata.images.values()))

    tiles = Tiles2D(0, len(image.coords["x"]), 0, len(image.coords["y"]), args.width)

    polygons = [
        poly
        for x_bounds, y_bounds in tqdm(tiles)
        for poly in run_patch(image, x_bounds, y_bounds, args.dapi, cellpose_patch())
    ]
    polygons = solve_conflicts(polygons)

    geo_df = gpd.GeoDataFrame({"geometry": polygons})
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
