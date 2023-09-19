import argparse

import geopandas as gpd
import numpy as np
import pandas as pd
import spatialdata
import xarray as xr
from anndata import AnnData
from shapely import affinity
from shapely.geometry import Polygon, box
from spatialdata.models import ShapesModel, TableModel
from tqdm import tqdm

from .._constants import ROI
from ..utils.tiling import Tiles2D
from ..utils.utils import _get_spatial_image
from .cellpose import cellpose_patch
from .utils import average, extract_polygons, solve_conflicts, to_chunk_mask

N_POSSIBLE_CHANNELS = [1, 2]


def run_patch(
    xarr: xr.DataArray,
    poly_ROI: Polygon | None,
    x_bounds: list[int],
    y_bounds: list[int],
    channels: list[int],
    method: callable,
    expand_radius: int = 0,
) -> np.ndarray:
    patch = xarr.sel(
        c=channels,
        x=slice(*x_bounds),
        y=slice(*y_bounds),
    ).values

    if poly_ROI is not None:
        patch_box = box(x_bounds[0], y_bounds[0], x_bounds[1], y_bounds[1])

        if not poly_ROI.intersects(patch_box):
            return []

        if not poly_ROI.contains(patch_box):
            patch = patch * to_chunk_mask(
                poly_ROI, x_bounds[0], y_bounds[0], x_bounds[1], y_bounds[1]
            )

    polygons = extract_polygons(method(patch), expand_radius)

    return [affinity.translate(p, x_bounds[0], y_bounds[0]) for p in polygons]


def main(args):
    sdata = spatialdata.read_zarr(args.path)

    image_key, image = _get_spatial_image(sdata)

    poly_ROI = sdata.shapes.get(ROI.KEY).geometry[0]

    tiles = Tiles2D(0, len(image.coords["x"]), 0, len(image.coords["y"]), args.width)

    assert (
        len(args.channels) in N_POSSIBLE_CHANNELS
    ), f"Provide one or two channel names. Found {len(args.channels)}"
    channels = [0, 0] if len(args.channels) == 1 else [1, 2]

    polygons = [
        poly
        for x_bounds, y_bounds in tqdm(tiles)
        for poly in run_patch(
            image,
            poly_ROI,
            x_bounds,
            y_bounds,
            args.channels,
            cellpose_patch(args.diameter, channels),
            args.expand_radius,
        )
    ]
    polygons = solve_conflicts(polygons)

    geo_df = gpd.GeoDataFrame(
        {
            "geometry": polygons,
            "x": [poly.centroid.x for poly in polygons],
            "y": [poly.centroid.y for poly in polygons],
        }
    )
    geo_df.index = image_key + geo_df.index.astype(str)

    geo_df = ShapesModel.parse(geo_df, transformations=image.transform)
    sdata.add_shapes("polygons", geo_df)

    mean_intensities = average(image, polygons)
    adata = AnnData(
        mean_intensities,
        dtype=mean_intensities.dtype,
        var=pd.DataFrame(index=image.c),
        obs=pd.DataFrame(index=geo_df.index),
    )

    adata.obsm["spatial"] = geo_df[["x", "y"]].values
    adata.obs["region"] = pd.Series("polygons", index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(image_key, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(image_key, index=adata.obs_names, dtype="category")
    adata.obs["cell_id"] = geo_df.index

    adata = TableModel.parse(
        adata,
        region_key="region",
        region="polygons",
        instance_key="cell_id",
    )

    sdata.table = adata


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
        "-c",
        "--channels",
        type=str,
        required=True,
        nargs="+",
        help="One or two channel names to be used for segmentation. If two channels, provide first the cytoplasm channel, and then the nucleus channel",
    )
    parser.add_argument(
        "-d",
        "--diameter",
        type=float,
        required=True,
        help="Expected cell diameter",
    )
    parser.add_argument(
        "-w",
        "--width",
        type=int,
        default=5000,
        help="Tile width",
    )
    parser.add_argument(
        "-e",
        "--expand_radius",
        type=int,
        default=0,
        help="Expand cell polygons by the provided radius",
    )

    main(parser.parse_args())
