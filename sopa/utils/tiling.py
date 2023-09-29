import logging
from math import ceil
from pathlib import Path

import dask.dataframe as dd
import geopandas as gpd
import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage
from shapely.geometry import Polygon, box
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import ShapesModel

from .._constants import ROI
from .utils import _get_item, to_intrinsic

log = logging.getLogger(__name__)


class Tiles1D:
    def __init__(self, xmin, xmax, tile_width, tile_overlap, int_coords):
        self.xmin, self.xmax = xmin, xmax
        self.delta = self.xmax - self.xmin

        self.tile_width = tile_width
        self.tile_overlap = tile_overlap
        self.int_coords = int_coords

        self._count = self.count()

    def count(self):
        if self.tile_width >= self.delta:
            return 1
        return ceil((self.delta - self.tile_overlap) / (self.tile_width - self.tile_overlap))

    def update(self, tile_width):
        self.tile_width = tile_width
        assert self._count == self.count()

    def tight_width(self):
        return ceil((self.delta + (self._count - 1) * self.tile_overlap) / self._count)

    def __getitem__(self, i):
        start_delta = i * (self.tile_width - self.tile_overlap)
        x0, x1 = self.xmin + start_delta, self.xmin + start_delta + self.tile_width

        return [int(x0), int(x1)] if self.int_coords else [x0, x1]


class Tiles2D:
    def __init__(
        self,
        sdata: SpatialData,
        element_name: str,
        tile_width: float | int,
        tile_overlap: float | int = 50,
    ):
        self.sdata = sdata
        self.element = sdata[element_name]

        if isinstance(self.element, MultiscaleSpatialImage):
            self.element = self.element["scale0"][element_name]

        if isinstance(self.element, SpatialImage) or isinstance(self.element, xr.DataArray):
            xmin, ymin = 0, 0
            xmax, ymax = len(self.element.coords["x"]), len(self.element.coords["y"])
            tight, int_coords = False, True
        elif isinstance(self.element, dd.DataFrame):
            xmin, ymin = self.element.x.min().compute(), self.element.y.min().compute()
            xmax, ymax = self.element.x.max().compute(), self.element.y.max().compute()
            tight, int_coords = True, False
        else:
            raise ValueError(f"Invalid element type: {type(self.element)}")

        self.tile_x = Tiles1D(xmin, xmax, tile_width, tile_overlap, int_coords)
        self.tile_y = Tiles1D(ymin, ymax, tile_width, tile_overlap, int_coords)

        self.tile_width = tile_width
        self.tile_overlap = tile_overlap
        self.tight = tight
        self.int_coords = int_coords
        self.roi = sdata.shapes[ROI.KEY] if ROI.KEY in sdata.shapes else None  # TODO: utils
        if self.roi is not None:
            self.roi = to_intrinsic(sdata, self.roi, element_name).geometry[0]

        assert self.tile_width > self.tile_overlap

        width_x = self.tile_x.tight_width()
        width_y = self.tile_y.tight_width()

        if self.tight:
            self.tile_width = max(width_x, width_y)
            self.tile_x.update(self.tile_width)
            self.tile_y.update(self.tile_width)

        self._ilocs = []

        for i in range(self.tile_x._count * self.tile_y._count):
            ix, iy = self.pair_indices(i)
            bounds = self.iloc(ix, iy)
            patch = box(*bounds)
            if self.roi is None or self.roi.intersects(patch):
                self._ilocs.append((ix, iy))

    def pair_indices(self, i: int) -> tuple[int, int]:
        iy, ix = divmod(i, self.tile_x._count)
        return ix, iy

    def iloc(self, ix: int, iy: int):
        xmin, xmax = self.tile_x[ix]
        ymin, ymax = self.tile_y[iy]
        return [xmin, ymin, xmax, ymax]

    def __getitem__(self, i) -> tuple[int, int, int, int]:
        """One tile a tuple representing its bounding box: (xmin, ymin, xmax, ymax)"""
        if isinstance(i, slice):
            start, stop, step = i.indices(len(self))
            return [self[i] for i in range(start, stop, step)]

        return self.iloc(*self._ilocs[i])

    def __len__(self):
        return len(self._ilocs)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def polygon(self, i: int) -> Polygon:
        square = box(*self[i])
        return square if self.roi is None else square.intersection(self.roi)

    @property
    def polygons(self) -> list[Polygon]:
        return [self.polygon(i) for i in range(len(self))]

    def write(self, overwrite: bool = True):
        geo_df = gpd.GeoDataFrame(
            {"geometry": self.polygons, "bounds": [self[i] for i in range(len(self))]}
        )
        geo_df = ShapesModel.parse(geo_df, transformations=self.element.attrs["transform"])
        self.sdata.add_shapes("patches", geo_df, overwrite=overwrite)

    def patchify_transcripts(
        self, baysor_dir: str, cell_key: str = None, unassigned_value: int | str = None
    ):
        import shapely
        from shapely.geometry import Point
        from tqdm import tqdm

        df = self.element

        if cell_key is not None and unassigned_value is not None:
            df[cell_key] = df[cell_key].replace(unassigned_value, 0)

        baysor_dir = Path(baysor_dir)

        log.info(f"Making {len(self)} sub-csv for Baysor")
        for i, polygon in enumerate(tqdm(self.polygons)):
            patch_dir = (baysor_dir / str(i)).absolute()
            patch_dir.mkdir(parents=True, exist_ok=True)
            patch_path = patch_dir / "transcripts.csv"

            tx0, ty0, tx1, ty1 = polygon.bounds
            where = (df.x >= tx0) & (df.x <= tx1) & (df.y >= ty0) & (df.y <= ty1)

            if polygon.area < box(*polygon.bounds).area:
                sub_df = df[where].compute()

                points = [Point(row) for row in sub_df[["x", "y"]].values]
                tree = shapely.STRtree(points)
                indices = tree.query(polygon, predicate="intersects")

                sub_df.iloc[indices].to_csv(patch_path)
            else:
                df[where].to_csv(patch_path, single_file=True)

        log.info(f"Patches saved in directory {baysor_dir}")
