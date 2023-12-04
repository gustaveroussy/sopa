import logging
from functools import partial
from math import ceil
from pathlib import Path

import dask.dataframe as dd
import geopandas as gpd
import pandas as pd
from dask.diagnostics import ProgressBar
from multiscale_spatial_image import MultiscaleSpatialImage
from shapely.geometry import Polygon, box
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from .._constants import ROI, SopaFiles, SopaKeys
from .._sdata import get_boundaries, get_item, get_spatial_image, to_intrinsic

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
    """Perform patching with overlap"""

    def __init__(
        self,
        sdata: SpatialData,
        element_name: str,
        patch_width: float | int,
        patch_overlap: float | int = 50,
    ):
        """
        Args:
            sdata: A `SpatialData` object
            element_name: Name of the element on with patches will be made
            patch_width: Width of the patches (in the unit of the coordinate system of the element)
            patch_overlap: Overlap width between the patches
        """
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
        """Index localization of one patch

        Args:
            i: The patch index

        Returns:
            A tuple `(index_x, index_y)` representing the 2D localization of the patch
        """
        iy, ix = divmod(i, self.patch_x._count)
        return ix, iy

    def iloc(self, ix: int, iy: int) -> list[int]:
        """Coordinates of the rectangle bounding box of the patch at the given indices

        Args:
            ix: Patch index in the x-axis
            iy: Patch indes in the y-axis

        Returns:
            A list `[xmin, ymin, xmax, ymax]` representing the bounding box of the patch
        """
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
        """Number of patches"""
        return len(self._ilocs)

    def __iter__(self):
        """Iterate over all patches (see `__getitem__`)"""
        for i in range(len(self)):
            yield self[i]

    def polygon(self, i: int) -> Polygon:
        """One patch as a shapely polygon. The polygon may not be just a square, if a ROI has been previously selected.

        Args:
            i: Patch index

        Returns:
            Polygon representing the patch
        """
        rectangle = box(*self[i])
        return rectangle if self.roi is None else rectangle.intersection(self.roi)

    @property
    def polygons(self) -> list[Polygon]:
        """All the patches as polygons

        Returns:
            List of `shapely` polygons
        """
        return [self.polygon(i) for i in range(len(self))]

    def write(self, overwrite: bool = True):
        """Save patches in `sdata.shapes["sopa_patches"]`

        Args:
            overwrite: Whether to overwrite patches if existing
        """
        geo_df = gpd.GeoDataFrame(
            {"geometry": self.polygons, SopaKeys.BOUNDS: [self[i] for i in range(len(self))]}
        )
        geo_df = ShapesModel.parse(
            geo_df, transformations=get_transformation(self.element, get_all=True)
        )
        self.sdata.add_shapes(SopaKeys.PATCHES, geo_df, overwrite=overwrite)

        log.info(f"{len(geo_df)} patches were saved in sdata['{SopaKeys.PATCHES}']")

    def patchify_transcripts(
        self,
        baysor_temp_dir: str,
        cell_key: str = None,
        unassigned_value: int | str = None,
        use_prior: bool = False,
    ) -> list[int]:
        """Patchification of the transcripts

        Args:
            baysor_temp_dir: Temporary directory where each patch will be stored. Note that each patch will have its own subdirectory.
            cell_key: Optional key of the transcript dataframe containing the cell IDs. This is useful if a prior segmentation has been run, assigning each transcript to a cell.
            unassigned_value: If `cell_key` has been provided, this corresponds to the value given in the 'cell ID' column for transcript that are not inside any cell.
            use_prior: Whether to use Cellpose as a prior segmentation for Baysor. If `True`, make sure you have already run Cellpose with Sopa, and no need to provide `cell_key` and `unassigned_value`. Note that, if you have MERFISH data, the prior has already been run, so just use `cell_key` and `unassigned_value`.

        Returns:
            A list of patches indices. Each index correspond to the name of a subdirectory inside `baysor_temp_dir`
        """
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

        if cell_key is None:
            cell_key = SopaKeys.BAYSOR_DEFAULT_CELL_KEY

        if unassigned_value is not None and unassigned_value != 0:
            self.df[cell_key] = self.df[cell_key].replace(unassigned_value, 0)

        if use_prior:
            prior_boundaries = self.sdata[SopaKeys.CELLPOSE_BOUNDARIES]
            _map_transcript_to_cell(self.sdata, cell_key, self.df, prior_boundaries)

        self._clean_directory()

        gdf = gpd.GeoDataFrame(geometry=self.patches_2d.polygons)
        with ProgressBar():
            self.df.map_partitions(partial(self._query_points_partition, gdf), meta=()).compute()

        log.info(f"Patches saved in directory {baysor_temp_dir}")
        return list(self.valid_indices())

    def _patch_path(self, index: int) -> Path:
        return self.baysor_temp_dir / str(index) / SopaFiles.BAYSOR_TRANSCRIPTS

    def _clean_directory(self):
        for index in range(len(self.patches_2d)):
            patch_path = self._patch_path(index)
            patch_path.parent.mkdir(parents=True, exist_ok=True)
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

    def _query_points_partition(self, gdf: gpd.GeoDataFrame, df: pd.DataFrame) -> pd.DataFrame:
        points_gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["x"], df["y"]))
        df_merged: gpd.GeoDataFrame = points_gdf.sjoin(gdf)

        for index, patch_df in df_merged.groupby("index_right"):
            patch_path = self._patch_path(index)
            patch_path.parent.mkdir(parents=True, exist_ok=True)
            patch_df = patch_df.drop(columns=["index_right", "geometry"])
            patch_df.to_csv(patch_path, mode="a", header=False, index=False)

    def _check_min_lines(self, path: str, n: int) -> bool:
        with open(path, "r") as f:
            for count, _ in enumerate(f):
                if count + 1 >= n:
                    return True
            return False


def _get_cell_id(gdf: gpd.GeoDataFrame, partition: pd.DataFrame, na_cells: int = 0) -> pd.Series:
    points_gdf = gpd.GeoDataFrame(
        partition, geometry=gpd.points_from_xy(partition["x"], partition["y"])
    )
    spatial_join = points_gdf.sjoin(gdf, how="left")
    spatial_join = spatial_join[~spatial_join.index.duplicated(keep="first")]
    cell_ids = (spatial_join["index_right"].fillna(-1) + 1 + na_cells).astype(int)

    return cell_ids


def _map_transcript_to_cell(
    sdata: SpatialData,
    cell_key: str,
    df: dd.DataFrame | None = None,
    geo_df: gpd.GeoDataFrame | None = None,
):
    if df is None:
        df = get_item(sdata, "points")

    if geo_df is None:
        geo_df = get_boundaries(sdata)

    geo_df = to_intrinsic(sdata, geo_df, df)
    geo_df = geo_df.reset_index()

    get_cell_id = partial(_get_cell_id, geo_df)

    if isinstance(df, dd.DataFrame):
        df[cell_key] = df.map_partitions(get_cell_id)
    else:
        raise ValueError(f"Invalid dataframe type: {type(df)}")
