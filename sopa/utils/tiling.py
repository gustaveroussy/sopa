from math import ceil

import numpy as np
import xarray as xr
from shapely.geometry import Polygon, box


class Tiles1D:
    def __init__(self, xmin, xmax, tile_width, tile_overlap, int_coords):
        self.xmin, self.xmax = xmin, xmax
        self.delta = self.xmax - self.xmin

        self.tile_width = tile_width
        self.tile_overlap = tile_overlap
        self.int_coords = int_coords

        self.count = self._count()

    def _count(self):
        if self.tile_width >= self.delta:
            return 1
        return ceil((self.delta - self.tile_overlap) / (self.tile_width - self.tile_overlap))

    def update(self, tile_width):
        self.tile_width = tile_width
        assert self.count == self._count()

    def tight_width(self):
        return ceil((self.delta + (self.count - 1) * self.tile_overlap) / self.count)

    def __getitem__(self, i):
        start_delta = i * (self.tile_width - self.tile_overlap)
        x0, x1 = self.xmin + start_delta, self.xmin + start_delta + self.tile_width

        return [int(x0), int(x1)] if self.int_coords else [x0, x1]


class Tiles2D:
    def __init__(
        self,
        tile_width: float | int,
        xmax: float | int,
        ymax: float | int,
        ymin: float | int = 0,
        xmin: float | int = 0,
        tile_overlap: float | int = 50,
        tight: bool = False,
        int_coords: bool = False,
    ):
        self.tile_x = Tiles1D(xmin, xmax, tile_width, tile_overlap, int_coords)
        self.tile_y = Tiles1D(ymin, ymax, tile_width, tile_overlap, int_coords)

        self.tile_width = tile_width
        self.tile_overlap = tile_overlap
        self.tight = tight
        self.int_coords = int_coords

        assert self.tile_width > self.tile_overlap

        self.compute()

    @classmethod
    def from_image(cls, image: xr.DataArray, tile_width: int, tile_overlap: int):
        xmax = len(image.coords["x"])
        ymax = len(image.coords["y"])
        return Tiles2D(tile_width, xmax, ymax, tile_overlap=tile_overlap, int_coords=True)

    def compute(self):
        width_x = self.tile_x.tight_width()
        width_y = self.tile_y.tight_width()

        self.count = self.tile_x.count * self.tile_y.count

        if self.tight:
            self.tile_width = max(width_x, width_y)
            self.tile_x.update(self.tile_width)
            self.tile_y.update(self.tile_width)

    def pair_indices(self, i: int) -> tuple[int, int]:
        iy, ix = divmod(i, self.tile_x.count)
        return ix, iy

    def __getitem__(self, i):
        if isinstance(i, slice):
            start, stop, step = i.indices(len(self))
            return [self[i] for i in range(start, stop, step)]
        ix, iy = self.pair_indices(i)

        xmin, xmax = self.tile_x[ix]
        ymin, ymax = self.tile_y[iy]

        return [xmin, ymin, xmax, ymax]

    def __len__(self):
        return self.count

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def coords_to_indices(self, coords: np.ndarray):
        assert self.tile_overlap == 0, "coords to indices with tile overlap is ambiguous"
        min_coords = np.array([self.tile_x.xmin, self.tile_y.xmin])

        return np.floor((coords - min_coords) / self.tile_width).clip(0).astype(int)

    def polygon(self, bounds: list[int]) -> Polygon:
        return box(*bounds)

    def polygons(self) -> list[Polygon]:
        return [self.polygon(self[i]) for i in range(len(self))]
