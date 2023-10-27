import logging
from functools import partial
from math import ceil
from pathlib import Path

import dask.dataframe as dd
import geopandas as gpd
from multiscale_spatial_image import MultiscaleSpatialImage
from shapely.geometry import Polygon, box
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from .._constants import ROI, SopaFiles, SopaKeys
from .._sdata import get_spatial_image, to_intrinsic
from .aggregate import map_transcript_to_cell
from .shapes import where_transcripts_inside_patch

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
        from tqdm import tqdm

        baysor_temp_dir = Path(baysor_temp_dir)
        df: dd.DataFrame = self.element

        if cell_key is not None and unassigned_value is not None:
            df[cell_key] = df[cell_key].replace(unassigned_value, 0)

        if cell_key is None:
            cell_key = SopaKeys.BAYSOR_DEFAULT_CELL_KEY

        prior_boundaries = None
        if use_prior:
            prior_boundaries = self.sdata[SopaKeys.CELLPOSE_BOUNDARIES]

        valid_indices = []

        for i, patch in enumerate(tqdm(self.polygons, desc="Splitting CSVs for Baysor")):
            patch_path: Path = baysor_temp_dir / str(i) / SopaFiles.BAYSOR_TRANSCRIPTS
            patch_path.parent.mkdir(parents=True, exist_ok=True)

            tx0, ty0, tx1, ty1 = patch.bounds
            sub_df = df[(df.x >= tx0) & (df.x <= tx1) & (df.y >= ty0) & (df.y <= ty1)]

            if patch.area < box(*patch.bounds).area:
                where_inside_patch = partial(where_transcripts_inside_patch, patch)
                sub_df = sub_df[sub_df.map_partitions(where_inside_patch)]

            if prior_boundaries is not None:
                map_transcript_to_cell(self.sdata, cell_key, sub_df, prior_boundaries)

            print(f"Computing CSV {i}")
            sub_df.compute().to_csv(patch_path)  # , single_file=True)
            print(f"Done {i}")

            if _check_min_lines(patch_path, 1000):
                valid_indices.append(i)
            else:
                log.info(f"Patch {i} has too few transcripts. Baysor will not be run on it.")

        log.info(f"Patches saved in directory {baysor_temp_dir}")
        return valid_indices


def _check_min_lines(path: str, n: int) -> bool:
    with open(path, "r") as f:
        for count, _ in enumerate(f):
            if count + 1 >= n:
                return True
        return False
