from math import ceil
from typing import Union

import numpy as np
from shapely.geometry import Polygon


class TileMaker1D:
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
        return ceil(
            (self.delta - self.tile_overlap) / (self.tile_width - self.tile_overlap)
        )

    def update(self, tile_width):
        self.tile_width = tile_width
        assert self.count == self._count()

    def tight_width(self):
        return ceil((self.delta + (self.count - 1) * self.tile_overlap) / self.count)

    def __getitem__(self, i):
        start_delta = i * (self.tile_width - self.tile_overlap)
        x0, x1 = self.xmin + start_delta, self.xmin + start_delta + self.tile_width

        return [int(x0), int(x1)] if self.int_coords else [x0, x1]


class TileMaker2D:
    def __init__(
        self,
        xmin: Union[float, int],
        xmax: Union[float, int],
        ymin: Union[float, int],
        ymax: Union[float, int],
        tile_width: Union[float, int],
        tile_overlap: Union[float, int] = 50,
        tight: bool = False,
        int_coords: bool = False,
    ):
        self.tile_x = TileMaker1D(xmin, xmax, tile_width, tile_overlap, int_coords)
        self.tile_y = TileMaker1D(ymin, ymax, tile_width, tile_overlap, int_coords)

        self.tile_width = tile_width
        self.tile_overlap = tile_overlap
        self.tight = tight
        self.int_coords = int_coords

        assert self.tile_width > self.tile_overlap

        self.compute()

    def compute(self):
        width_x = self.tile_x.tight_width()
        width_y = self.tile_y.tight_width()

        self.count = self.tile_x.count * self.tile_y.count

        if self.tight:
            self.tile_width = max(width_x, width_y)
            self.tile_x.update(self.tile_width)
            self.tile_y.update(self.tile_width)

    def __getitem__(self, i):
        iy, ix = divmod(i, self.tile_x.count)
        return [self.tile_x[ix], self.tile_y[iy]]

    def __len__(self):
        return self.count

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def polygon(self, i):
        (x0, x1), (y0, y1) = self[i]
        return Polygon([(x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0)])

    def polygons(self):
        return [self.polygon(i) for i in range(len(self))]

    def coords_to_indices(self, coords: np.ndarray):
        assert (
            self.tile_overlap == 0
        ), "coords to indices with tile overlap is ambiguous"
        min_coords = np.array([self.tile_x.xmin, self.tile_y.xmin])

        return np.floor((coords - min_coords) / self.tile_width).clip(0).astype(int)
