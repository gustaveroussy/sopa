from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity

from ...segmentation.nucleus import pad
from . import write_experiment, write_np_ome_tif, write_polygons, write_transcripts

path = Path("/Volumes/T7_Quentin/data/vizgen/example/region_0")

res_dir = path / "explorer"
res_dir.mkdir(exist_ok=True)

### JSON

write_experiment(res_dir / "experiment.xenium", "A0", "Breast", "Breast")

### IMAGE

# image = path / "images" / "mosaic_DAPI_z3.tif"
# image = np.zeros((28250, 45030, 1), dtype="uint8")

# write_np_ome_tif(
#     res_dir / "morphology.ome.tif", image, channel_names=["DAPI"], subresolutions=7
# )

### POLY

microns_to_pixels = np.loadtxt(path / "images" / "micron_to_mosaic_pixel_transform.csv")
# matrix = np.concatenate(
#     [microns_to_pixels[0, :2], microns_to_pixels[1, :2], microns_to_pixels[:2, 2]]
# )

# gdf = gpd.read_parquet(path / "cell_boundaries.parquet")
# polygons = gdf.Geometry.map(lambda mp: mp.geoms[0])
# polygons = [affinity.affine_transform(p, matrix) for p in polygons]

# coordinates = np.stack([pad(p, 3, 13) for p in polygons])
# coordinates /= 4.705882

# write_polygons(res_dir / "cells.zarr.zip", coordinates)


# ### TRANSCRIPTS

df = pd.read_csv(path / "detected_transcripts.csv")
x = "global_x"
y = "global_y"
xy = np.concatenate([df[[x, y]].values, np.ones((len(df), 1))], axis=1)
df[[x, y]] = (microns_to_pixels @ xy.T / 4.705882).T[:, :2]

print(df[[x, y]].min())

write_transcripts(res_dir / "transcripts.zarr.zip", df, x="global_x", y="global_y")
