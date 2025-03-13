import logging
from math import ceil

import dask.dataframe as dd
import geopandas as gpd
import numpy as np
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, box
from shapely.ops import unary_union
from spatialdata import SpatialData
from spatialdata.models import ShapesModel, SpatialElement
from spatialdata.transformations import get_transformation
from xarray import DataArray, DataTree

from .._constants import SopaKeys
from ..utils import add_spatial_element, to_intrinsic

log = logging.getLogger(__name__)


class Patches1D:
    def __init__(
        self,
        xmin: float | int,
        xmax: float | int,
        patch_width: float | int,
        patch_overlap: float | int,
        tight: bool,
        int_coords: bool,
    ):
        self.xmin, self.xmax = xmin, xmax
        self.delta = self.xmax - self.xmin

        self.patch_width = patch_width
        self.patch_overlap = patch_overlap
        self.int_coords = int_coords

        self._count = self.count()
        if tight:
            self.patch_width = self.tight_width()
            assert (
                self._count == self.count()
            ), f"Invalid patching with {self.delta=}, {self.patch_width=} and {self.patch_overlap=}"

    def count(self):
        if self.patch_width >= self.delta:
            return 1
        return ceil((self.delta - self.patch_overlap) / (self.patch_width - self.patch_overlap))

    def update(self, patch_width):
        self.patch_width = patch_width
        assert self._count == self.count()

    def tight_width(self):
        width = (self.delta + (self._count - 1) * self.patch_overlap) / self._count
        return ceil(width) if self.int_coords else width + 1e-4  # add an epsilon to avoid floating point errors (#214)

    def __getitem__(self, i):
        start_delta = i * (self.patch_width - self.patch_overlap)
        x0, x1 = self.xmin + start_delta, self.xmin + start_delta + self.patch_width

        return [int(x0), int(x1)] if self.int_coords else [x0, x1]


class Patches2D:
    """
    Compute 2D-patches with overlaps. This can be done on an image or a DataFrame.

    Attributes:
        polygons (list[Polygon]): List of `shapely` polygons representing the patches
        bboxes (np.ndarray): Array of shape `(n_patches, 4)` containing the (xmin, ymin, xmax, ymax) coordinates of the patches bounding boxes
        ilocs (np.ndarray): Array of shape `(n_patches, 2)` containing the (x,y) indices of the patches
    """

    polygons: list[Polygon]
    bboxes: np.ndarray
    ilocs: np.ndarray

    def __init__(
        self,
        sdata: SpatialData,
        element: SpatialElement | str,
        patch_width: float | int | None,
        patch_overlap: float | int = 50,
    ):
        """
        Args:
            sdata: A `SpatialData` object
            element: SpatialElement or name of the element on with patches will be made (image or points)
            patch_width: Width of the patches (in the unit of the coordinate system of the element). If `None`, only one patch will be created
            patch_overlap: Overlap width between the patches
        """
        patch_width = float("inf") if (patch_width is None or patch_width == -1) else patch_width

        assert patch_width > patch_overlap, f"Argument {patch_width=} must be greater than {patch_overlap=}"

        self.sdata = sdata
        self.element = sdata[element] if isinstance(element, str) else element
        self.original_element = self.element  # in case the element is a DataTree

        if isinstance(self.element, DataTree):
            self.element = next(iter(self.element["scale0"].values()))

        if isinstance(self.element, DataArray):
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

        self.roi = None
        if SopaKeys.ROI in sdata.shapes:
            geo_df = to_intrinsic(sdata, sdata[SopaKeys.ROI], self.original_element)

            assert all(
                isinstance(geom, Polygon) for geom in geo_df.geometry
            ), f"All sdata['{SopaKeys.ROI}'] geometries must be polygons"

            self.roi = unary_union(geo_df.geometry)  # merge polygons into one multi-polygon
            assert isinstance(self.roi, (Polygon, MultiPolygon)), f"Invalid ROI type: {type(self.roi)}"

        self._init_patches()

    def _init_patches(self):
        self.ilocs, self.polygons, self.bboxes = [], [], []

        for i in range(self.patch_x._count * self.patch_y._count):
            self._try_add_patch(i)

        self.ilocs = np.array(self.ilocs)
        self.bboxes = np.array(self.bboxes)

    def _try_add_patch(self, i: int):
        """Check that the patch is valid, and, if valid, add it to the list of valid patches"""
        iy, ix = divmod(i, self.patch_x._count)
        bounds = self._bbox_iloc(ix, iy)
        patch = box(*bounds)

        if self.roi is not None and not self.roi.intersects(patch):
            return

        patch = patch if self.roi is None else patch.intersection(self.roi)

        if isinstance(patch, GeometryCollection):
            geoms = [geom for geom in patch.geoms if isinstance(geom, Polygon)]
            if not geoms:
                return
            patch = max(geoms, key=lambda polygon: polygon.area)

        if not isinstance(patch, (Polygon, MultiPolygon)):
            return

        self.polygons.append(patch)
        self.ilocs.append((ix, iy))
        self.bboxes.append(bounds)

    def __repr__(self):
        return f"{self.__class__.__name__} object with {len(self)} patches"

    @property
    def shape(self) -> tuple[int, int]:
        return (self.patch_y._count, self.patch_x._count)

    def _bbox_iloc(self, ix: int, iy: int) -> list[int]:
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

    def __len__(self):
        """Number of patches"""
        return len(self.bboxes)

    def write(self, *args, **kwargs):
        log.warning(
            "Patches2D.write is deprecated and will be removed in sopa==2.1.0. Use Patches2D.add_shapes instead"
        )
        self.add_shapes(*args, **kwargs)

    def add_shapes(self, key_added: str | None = None, overwrite: bool = True) -> gpd.GeoDataFrame:
        key_added = key_added or SopaKeys.PATCHES
        geo_df = self.as_geodataframe()

        add_spatial_element(self.sdata, key_added, geo_df, overwrite=overwrite)

        log.info(f"Added {len(geo_df)} patche(s) to sdata['{key_added}']")

        return geo_df

    def as_geodataframe(self) -> gpd.GeoDataFrame:
        geo_df = gpd.GeoDataFrame(
            {
                "geometry": self.polygons,
                SopaKeys.BOUNDS: self.bboxes.tolist(),
                SopaKeys.PATCHES_ILOCS: self.ilocs.tolist(),
            }
        )

        return ShapesModel.parse(geo_df, transformations=get_transformation(self.element, get_all=True).copy())

    def patchify_transcripts(self, *args, **kwargs):
        raise NameError(
            "Patches2D.patchify_transcripts is deprecated and will be removed in sopa==2.1.0. Use `sopa.make_transcript_patches` instead"
        )

    def patchify_centroids(self, *args, **kwargs) -> list[int]:
        raise NameError(
            "Patches2D.patchify_centroids is deprecated and will be removed in sopa==2.1.0. Use `sopa.make_transcript_patches` instead"
        )
