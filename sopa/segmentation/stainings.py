import logging
from pathlib import Path

import numpy as np
import zarr
from shapely import affinity
from shapely.geometry import Polygon, box
from spatialdata import SpatialData
from tqdm import tqdm

from .._constants import SopaKeys
from .._sdata import get_spatial_image
from . import shapes
from .patching import Patches2D

log = logging.getLogger(__name__)


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

        self.image_key, self.image = get_spatial_image(sdata, return_key=True)

        assert np.isin(
            channels, self.image.c
        ).all(), f"Channel names must be a subset of: {', '.join(self.image.c)}"

    def _run_patch(
        self,
        patch: Polygon,
    ) -> list[Polygon]:
        bounds = [int(x) for x in patch.bounds]

        image = self.image.sel(
            c=self.channels,
            x=slice(bounds[0], bounds[2]),
            y=slice(bounds[1], bounds[3]),
        ).values

        if patch.area < box(*bounds).area:
            image = image * shapes.rasterize(patch, image.shape[1:], bounds)

        cells = shapes.geometrize(self.method(image))

        return [affinity.translate(cell, *bounds[:2]) for cell in cells]

    def write_patch_cells(self, patch_dir: str, patch_index: int):
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
        patch_width: int,
        patch_overlap: int,
    ) -> list[Polygon]:
        self.patches = Patches2D(self.sdata, self.image_key, patch_width, patch_overlap)

        cells = [
            cell
            for patch in tqdm(self.patches.polygons, desc="Run on patches")
            for cell in self._run_patch(patch)
        ]
        cells = shapes.solve_conflicts(cells)
        return cells

    @classmethod
    def read_patches_cells(cls, patch_dir: str) -> list[Polygon]:
        cells = []

        files = [f for f in Path(patch_dir).iterdir() if f.suffix == ".zip"]
        for file in tqdm(files, desc="Reading patches"):
            z = zarr.open(file, mode="r")
            for _, coords_zarr in z.arrays():
                cells.append(Polygon(coords_zarr[:]))

        log.info(f"Found {len(cells)} total cells")
        return cells
