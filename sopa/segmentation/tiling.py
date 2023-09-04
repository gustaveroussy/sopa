from math import ceil
from typing import Union


class TileMaker1D:
    def __init__(self, xmin, xmax, tile_width, tile_overlap, unit):
        self.xmin, self.xmax = xmin, xmax
        self.delta = self.xmax - self.xmin

        self.tile_width = tile_width
        self.tile_overlap = tile_overlap
        self.unit = unit

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

    def coordinates(self, i):
        start_delta = i * (self.tile_width - self.tile_overlap)
        x0, x1 = self.xmin + start_delta, self.xmin + start_delta + self.tile_width

        return [int(x0), int(x1)] if self.unit == "pixel" else [x0, x1]


class TileMaker:
    VALID_UNITS = ["pixel", "micron"]

    def __init__(
        self,
        xmin: Union[float, int],
        xmax: Union[float, int],
        ymin: Union[float, int],
        ymax: Union[float, int],
        tile_width: Union[float, int],
        tile_overlap: Union[float, int] = 50,
        unit: str = "pixel",
    ):
        self.tile_x = TileMaker1D(xmin, xmax, tile_width, tile_overlap, unit)
        self.tile_y = TileMaker1D(ymin, ymax, tile_width, tile_overlap, unit)

        self.tile_width = tile_width
        self.tile_overlap = tile_overlap
        self.unit = unit

        assert (
            self.unit in self.VALID_UNITS
        ), f"Argument {self.unit} is not a valid unit. It must be one of: {', '.join(self.VALID_UNITS)}"

        assert self.tile_width > self.tile_overlap

        self.compute()

    def compute(self):
        width_x = self.tile_x.tight_width()
        width_y = self.tile_y.tight_width()

        self.count = self.tile_x.count * self.tile_y.count

        self.tile_width = max(width_x, width_y)
        self.tile_x.update(self.tile_width)
        self.tile_y.update(self.tile_width)

    def coordinates(self, i):
        iy, ix = divmod(i, self.tile_x.count)
        return [self.tile_x.coordinates(ix), self.tile_y.coordinates(iy)]

    def __len__(self):
        return self.count

    def __iter__(self):
        for i in range(self.count):
            yield self.coordinates(i)
