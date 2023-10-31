import logging
from functools import partial
from math import ceil
from pathlib import Path

import dask.dataframe as dd
import geopandas as gpd
import pandas as pd
import shapely
from dask.diagnostics import ProgressBar
from multiscale_spatial_image import MultiscaleSpatialImage
from shapely.geometry import Point, Polygon, box
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from .._constants import ROI, SopaFiles, SopaKeys
from .._sdata import get_spatial_image, to_intrinsic
from .aggregate import map_transcript_to_cell

log = logging.getLogger(__name__)


class Patches1D:
    def __init__(self, xmin, xmax, patch_width, patch_overlap, tight, int_coords):
        self.xmin, self.xmax = xmin, xmax
        self.delta = self.xmax - self.xmin

        self.patch_width = patch_width
        self.patch_overlap = patch_overlap
        self.int_coords = int_coords

        self._count = self.count()
        if tight:
            self.patch_width = self.tight_width()
            assert self._count == self.count()

    def count(self):
        if self.patch_width >= self.delta:
            return 1
        return ceil((self.delta - self.patch_overlap) / (self.patch_width - self.patch_overlap))

    def update(self, patch_width):
        self.patch_width = patch_width
        assert self._count == self.count()

    def tight_width(self):
        width = (self.delta + (self._count - 1) * self.patch_overlap) / self._count
        return ceil(width) if self.int_coords else width

    def __getitem__(self, i):
        start_delta = i * (self.patch_width - self.patch_overlap)
        x0, x1 = self.xmin + start_delta, self.xmin + start_delta + self.patch_width

        return [int(x0), int(x1)] if self.int_coords else [x0, x1]


class Patches2D:
    def __init__(
        self,
        sdata: SpatialData,
        element_name: str,
        patch_width: float | int,
        patch_overlap: float | int = 50,
    ):
        self.sdata = sdata
        self.element_name = element_name
        self.element = sdata[element_name]

        if isinstance(self.element, MultiscaleSpatialImage):
            self.element = get_spatial_image(sdata, element_name)

        if isinstance(self.element, SpatialImage):
            xmin, ymin = 0, 0
            xmax, ymax = len(self.element.coords["x"]), len(self.element.coords["y"])
            tight, int_coords = False, True
        elif isinstance(self.element, dd.DataFrame):
            xmin, ymin = self.element.x.min().compute(), self.element.y.min().compute()
            xmax, ymax = self.element.x.max().compute(), self.element.y.max().compute()
            tight, int_coords = True, False
        else:
            raise ValueError(f"Invalid element type: {type(self.element)}")

        self.patch_x = Patches1D(xmin, xmax, patch_width, patch_overlap, tight, int_coords)
        self.patch_y = Patches1D(ymin, ymax, patch_width, patch_overlap, tight, int_coords)

        self.roi = sdata.shapes[ROI.KEY] if ROI.KEY in sdata.shapes else None
        if self.roi is not None:
            self.roi = to_intrinsic(sdata, self.roi, element_name).geometry[0]

        self._ilocs = []

        for i in range(self.patch_x._count * self.patch_y._count):
            ix, iy = self.pair_indices(i)
            bounds = self.iloc(ix, iy)
            patch = box(*bounds)
            if self.roi is None or self.roi.intersects(patch):
                self._ilocs.append((ix, iy))

    def pair_indices(self, i: int) -> tuple[int, int]:
        iy, ix = divmod(i, self.patch_x._count)
        return ix, iy

    def iloc(self, ix: int, iy: int):
        xmin, xmax = self.patch_x[ix]
        ymin, ymax = self.patch_y[iy]
        return [xmin, ymin, xmax, ymax]

    def __getitem__(self, i) -> tuple[int, int, int, int]:
        """One patche bounding box: (xmin, ymin, xmax, ymax)"""
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
        rectangle = box(*self[i])
        return rectangle if self.roi is None else rectangle.intersection(self.roi)

    @property
    def polygons(self) -> list[Polygon]:
        return [self.polygon(i) for i in range(len(self))]

    def write(self, overwrite: bool = True):
        geo_df = gpd.GeoDataFrame(
            {"geometry": self.polygons, SopaKeys.BOUNDS: [self[i] for i in range(len(self))]}
        )
        geo_df = ShapesModel.parse(
            geo_df, transformations=get_transformation(self.element, get_all=True)
        )
        self.sdata.add_shapes(SopaKeys.PATCHES, geo_df, overwrite=overwrite)

        log.info(f"{len(geo_df)} patches where saved in sdata['{SopaKeys.PATCHES}']")

    def patchify_transcripts(
        self,
        baysor_temp_dir: str,
        cell_key: str = None,
        unassigned_value: int | str = None,
        use_prior: bool = False,
    ) -> list[int]:
        return BaysorPatches(self, self.element).write(
            baysor_temp_dir, cell_key, unassigned_value, use_prior
        )


class BaysorPatches:
    def __init__(self, patches_2d: Patches2D, df: dd.DataFrame):
        self.patches_2d = patches_2d
        self.df = df
        self.sdata = self.patches_2d.sdata

    def write(
        self,
        baysor_temp_dir: str,
        cell_key: str = None,
        unassigned_value: int | str = None,
        use_prior: bool = False,
    ):
        log.info("Writing sub-CSV for baysor")
        self.baysor_temp_dir = Path(baysor_temp_dir)
        self._clean_directory()

        if cell_key is not None and unassigned_value is not None:
            self.df[cell_key] = self.df[cell_key].replace(unassigned_value, 0)

        if cell_key is None:
            cell_key = SopaKeys.BAYSOR_DEFAULT_CELL_KEY

        if use_prior:
            prior_boundaries = self.sdata[SopaKeys.CELLPOSE_BOUNDARIES]
            map_transcript_to_cell(self.sdata, cell_key, self.df, prior_boundaries)

        tree = shapely.STRtree(self.patches_2d.polygons)
        with ProgressBar():
            self.df.map_partitions(partial(self._query_points_partition, tree)).compute()

        log.info(f"Patches saved in directory {baysor_temp_dir}")
        return list(self.valid_indices())

    def _patch_path(self, index: int) -> Path:
        return self.baysor_temp_dir / str(index) / SopaFiles.BAYSOR_TRANSCRIPTS

    def _clean_directory(self):
        for index in range(len(self.patches_2d)):
            patch_path = self._patch_path(index)
            if patch_path.exists():
                patch_path.unlink()
            pd.DataFrame(columns=self.df.columns).to_csv(patch_path, index=False)

    def valid_indices(self):
        for index in range(len(self.patches_2d)):
            patch_path = self._patch_path(index)
            if self._check_min_lines(patch_path, 1000):
                yield index
            else:
                log.info(f"Patch {index} has < 1000 transcripts. Baysor will not be run on it.")

    def _query_points_partition(
        self, tree: shapely.STRtree, df: pd.DataFrame, partition_info=None
    ) -> pd.DataFrame:
        points = [Point(xy) for xy in zip(df.x, df.y)]
        df_query = pd.DataFrame(tree.query(points).T, columns=["point_index", "patch_index"])
        df_merged = df.merge(df_query, left_index=True, right_on="point_index", how="right")

        if partition_info is not None:
            for index, patch_df in df_merged.groupby("patch_index"):
                patch_path = self._patch_path(index)
                patch_path.parent.mkdir(parents=True, exist_ok=True)
                patch_df = patch_df.drop(columns=["patch_index", "point_index"])
                patch_df.to_csv(patch_path, mode="a", header=False, index=False)

    def _check_min_lines(self, path: str, n: int) -> bool:
        with open(path, "r") as f:
            for count, _ in enumerate(f):
                if count + 1 >= n:
                    return True
            return False
