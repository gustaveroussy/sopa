from pathlib import Path

import numpy as np
import zarr
from shapely import affinity
from shapely.geometry import Polygon, box
from spatialdata import SpatialData
from tqdm import tqdm

from .._constants import SopaKeys
from ..utils.tiling import Tiles2D
from ..utils.utils import _get_spatial_image
from . import shapes


class StainingSegmentation:
    def __init__(
        self,
        sdata: SpatialData,
        method: callable,
        channels: list[str],
    ):
        self.sdata = sdata
        self.method = method
        self.channels = channels

        self.image_key, self.image = _get_spatial_image(sdata)

        assert np.isin(
            channels, self.image.c
        ).all(), f"Channel names must be a subset of: {', '.join(self.image.c)}"

    def _run_patch(
        self,
        polygon: Polygon,
    ) -> list[Polygon]:
        bounds = [int(x) for x in polygon.bounds]

        patch = self.image.sel(
            c=self.channels,
            x=slice(bounds[0], bounds[2]),
            y=slice(bounds[1], bounds[3]),
        ).values

        if polygon.area < box(*bounds).area:
            bounds = shapes.update_bounds(bounds, patch.shape[1:])
            patch = patch * shapes.rasterize(polygon, bounds)

        polygons = shapes.geometrize(self.method(patch))

        return [affinity.translate(p, *bounds[:2]) for p in polygons]

    def write_patch_polygons(self, patch_dir: str, patch_index: int):
        patch = self.sdata[SopaKeys.PATCHES].geometry[patch_index]
        cells = self._run_patch(patch)

        patch_file = Path(patch_dir) / f"{patch_index}.zarr.zip"

        with zarr.ZipStore(patch_file, mode="w") as store:
            g = zarr.group(store=store)

            for i, cell in enumerate(cells):
                coords = np.array(cell.exterior.coords)
                g.array(f"cell_{i}", coords, dtype=coords.dtype, chunks=coords.shape)

    def run_patches(
        self,
        tile_width: int,
        tile_overlap: int,
    ) -> list[Polygon]:
        self.tiles = Tiles2D(self.sdata, self.image_key, tile_width, tile_overlap)

        cells = [cell for patch in tqdm(self.tiles.polygons) for cell in self._run_patch(patch)]
        cells = shapes.solve_conflicts(cells)
        return cells
