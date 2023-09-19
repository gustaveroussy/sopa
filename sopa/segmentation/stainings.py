from pathlib import Path

import numpy as np
import zarr
from shapely import affinity
from shapely.geometry import Polygon
from spatialdata import SpatialData
from tqdm import tqdm

from .._constants import ROI
from ..utils.tiling import Tiles2D
from ..utils.utils import _get_spatial_image
from .shapes import extract_polygons, solve_conflicts, to_chunk_mask


class StainingSegmentation:
    def __init__(
        self,
        sdata: SpatialData,
        method: callable,
        channels: list[str],
        tile_width: int,
        tile_overlap: int,
        expand_radius: int = 0,
    ):
        self.sdata = sdata
        self.method = method
        self.channels = channels
        self.expand_radius = expand_radius

        self.image_key, self.image = _get_spatial_image(sdata)

        self.tiles = Tiles2D.from_image(self.image, tile_width, tile_overlap)

        assert np.isin(
            channels, self.image.c
        ).all(), f"Channel names must be a subset of: {', '.join(self.image.c)}"

    @property
    def poly_ROI(self) -> Polygon | None:
        if ROI.KEY in self.sdata.shapes:
            return self.sdata.shapes[ROI.KEY].geometry[0]
        return None

    def _run_patch(
        self,
        bounds: list[int],
    ) -> list[Polygon]:
        patch = self.image.sel(
            c=self.channels,
            x=slice(bounds[0], bounds[2]),
            y=slice(bounds[1], bounds[3]),
        ).values

        if self.poly_ROI is not None:
            patch_box = self.tiles.polygon(bounds)

            if not self.poly_ROI.intersects(patch_box):
                return []

            if not self.poly_ROI.contains(patch_box):
                patch = patch * to_chunk_mask(self.poly_ROI, bounds)

        polygons = extract_polygons(self.method(patch), self.expand_radius)

        return [affinity.translate(p, *bounds[:2]) for p in polygons]

    def patch_polygons_paths(self, directory: str):
        directory = Path(directory)
        return [directory / f"{i}.zarr.zip" for i in range(len(self.tiles))]

    def write_patch_polygons(self, path: str, index: int):
        polygons = self._run_patch(self.tiles[index])

        with zarr.ZipStore(path, mode="w") as store:
            g = zarr.group(store=store)

            for i, polygon in enumerate(polygons):
                coords = np.array(polygon.exterior.coords)
                g.array(f"polygon_{i}", coords, dtype=coords.dtype, chunks=coords.shape)

    def aggregate_patch_polygons(self, directory: str) -> list[Polygon]:
        polygons = []

        for path in tqdm(self.patch_polygons_paths(directory)):
            z = zarr.open(path, mode="r")
            for _, coords_zarr in z.arrays():
                polygons.append(Polygon(coords_zarr[:]))

        polygons = solve_conflicts(polygons)
        return polygons

    def run_patches(self) -> list[Polygon]:
        polygons = [poly for bounds in tqdm(self.tiles) for poly in self._run_patch(bounds)]
        polygons = solve_conflicts(polygons)
        return polygons
