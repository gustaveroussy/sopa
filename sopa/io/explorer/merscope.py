import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from shapely import affinity
from shapely.geometry import shape

from ...segmentation.nucleus import pad
from . import (
    write_experiment,
    write_groups,
    write_multiscale,
    write_polygons,
    write_transcripts,
)

path = Path("/Volumes/T7_Quentin/data/vizgen/example/region_2")


res_dir = path / "explorer"
res_dir.mkdir(exist_ok=True)

image_path = res_dir / "morphology.ome.tif"
cells_path = res_dir / "cells.zarr.zip"
transcripts_path = res_dir / "transcripts.zarr.zip"
experiment_path = res_dir / "experiment.xenium"
analysis_path = res_dir / "analysis.zarr.zip"

### JSON

write_experiment(experiment_path, "B2", "Breast", "Breast")

### IMAGE

if not image_path.exists():
    from dask_image.imread import imread
    from multiscale_spatial_image import to_multiscale
    from spatial_image import to_spatial_image

    im = imread(path / "images" / "mosaic_DAPI_z3.tif")
    im = to_spatial_image(im, dims=["c", "y", "x"], c_coords=["DAPI"])
    multiscale = to_multiscale(im, [2, 2, 2, 2, 2])

    write_multiscale(image_path, multiscale)

### TABLE

adata = anndata.read_h5ad(path / "baysor" / "adata_annotated_filtered.h5ad")
valid_cell = adata.obs["cell_num"].values

if not analysis_path.exists():
    write_groups(analysis_path, adata.obs)

### POLY

microns_to_pixels = np.loadtxt(path / "images" / "micron_to_mosaic_pixel_transform.csv")

if not cells_path.exists():
    matrix = np.concatenate(
        [microns_to_pixels[0, :2], microns_to_pixels[1, :2], microns_to_pixels[:2, 2]]
    )

    # gdf = gpd.read_parquet(path / "cell_boundaries.parquet")
    # polygons = gdf.Geometry.map(lambda mp: mp.geoms[0])
    with open(path / "baysor" / "segmentation_polygons.json") as f:
        d = json.load(f)
        d = {c["cell"]: c for c in d["geometries"]}
    polygons = [shape(d[cell_id]) for cell_id in valid_cell]
    polygons = [affinity.affine_transform(p, matrix) for p in polygons]

    coordinates = np.stack([pad(p, 3, 13) for p in polygons])
    coordinates /= 4.705882

    write_polygons(cells_path, coordinates)


# ### TRANSCRIPTS

if not transcripts_path.exists():
    df = pd.read_csv(path / "detected_transcripts.csv")
    x = "global_x"
    y = "global_y"
    xy = np.concatenate([df[[x, y]].values, np.ones((len(df), 1))], axis=1)
    df[[x, y]] = (microns_to_pixels @ xy.T / 4.705882).T[:, :2]

    print(df[[x, y]].min())

    write_transcripts(transcripts_path, df, x="global_x", y="global_y")
