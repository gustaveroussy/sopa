import argparse

import geopandas as gpd
import pandas as pd
import spatialdata
from anndata import AnnData
from spatialdata.models import ShapesModel, TableModel

from .cellpose import cellpose_patch
from .shapes import average
from .stainings import StainingSegmentation


def main(args):
    sdata = spatialdata.read_zarr(args.sdata_path)

    method = cellpose_patch(args.diameter, args.channels)
    segmentation = StainingSegmentation(
        sdata, method, args.channels, args.tile_width, args.tile_overlap, args.expand_radius
    )
    polygons = segmentation.run_patches()

    geo_df = gpd.GeoDataFrame(
        {
            "geometry": polygons,
            "x": [poly.centroid.x for poly in polygons],
            "y": [poly.centroid.y for poly in polygons],
        }
    )
    geo_df.index = segmentation.image_key + geo_df.index.astype(str)

    geo_df = ShapesModel.parse(geo_df, transformations=segmentation.image.transform)
    sdata.add_shapes("polygons", geo_df)

    mean_intensities = average(segmentation.image, polygons)
    adata = AnnData(
        mean_intensities,
        dtype=mean_intensities.dtype,
        var=pd.DataFrame(index=segmentation.image.c),
        obs=pd.DataFrame(index=geo_df.index),
    )

    adata.obsm["spatial"] = geo_df[["x", "y"]].values
    adata.obs["region"] = pd.Series("polygons", index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(segmentation.image_key, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(
        segmentation.image_key, index=adata.obs_names, dtype="category"
    )
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
        "-s",
        "--sdata_path",
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
        "-tw",
        "--tile_width",
        type=int,
        default=5000,
        help="Tile width (pixels)",
    )
    parser.add_argument(
        "-to",
        "--tile_overlap",
        type=int,
        default=50,
        help="Tile overlap (pixels)",
    )
    parser.add_argument(
        "-e",
        "--expand_radius",
        type=int,
        default=0,
        help="Expand cell polygons by the provided radius",
    )

    main(parser.parse_args())
