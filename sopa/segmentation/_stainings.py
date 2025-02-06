import logging
from functools import partial
from pathlib import Path
from typing import Callable, Iterable

import geopandas as gpd
import numpy as np
import pandas as pd
import spatialdata
from scipy.ndimage import gaussian_filter
from shapely.geometry import Polygon, box
from skimage import exposure
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation

from .. import settings
from .._constants import SopaKeys
from ..utils import add_spatial_element, get_spatial_image
from . import shapes

log = logging.getLogger(__name__)


class StainingSegmentation:
    def __init__(
        self,
        sdata: SpatialData,
        method: Callable,
        channels: list[str] | str,
        image_key: str | None = None,
        min_area: float = 0,
        clip_limit: float = 0.2,
        clahe_kernel_size: int | Iterable[int] | None = None,
        gaussian_sigma: float = 1,
    ):
        """Generalized staining-based segmentation class

        Args:
            sdata: A `SpatialData` object
            method: A segmentation `callable` whose input is an image of shape `(C, Y, X)` and output is a cell mask of shape `(Y, X)`. Each mask value `>0` represent a unique cell ID. Such callables can be found in `sopa.segmentation.methods`.
            channels: One or a list of channel names used for segmentation. If only one channel is provided, the image given to the `method` will be of shape `(1, Y, X)`.
            image_key: Optional key of `sdata` containing the image (no needed if there is only one image)
            min_area: Minimum area (in pixels^2) for a cell to be kept
            clip_limit: Parameter for skimage.exposure.equalize_adapthist (applied before running cellpose)
            clahe_kernel_size: Parameter for skimage.exposure.equalize_adapthist (applied before running cellpose)
            gaussian_sigma: Parameter for scipy gaussian_filter (applied before running cellpose)
        """
        assert SopaKeys.PATCHES in sdata.shapes, "Run `sopa.make_image_patches` before running segmentation"

        self.patches_gdf: gpd.GeoDataFrame = sdata[SopaKeys.PATCHES]
        self.method = method
        self.channels = [channels] if isinstance(channels, str) else channels

        self.min_area = min_area
        self.clip_limit = clip_limit
        self.clahe_kernel_size = clahe_kernel_size
        self.gaussian_sigma = gaussian_sigma

        self.image_key, self.image = get_spatial_image(sdata, key=image_key, return_key=True)

        image_channels = self.image.coords["c"].values
        assert np.isin(
            channels, image_channels
        ).all(), f"Channel names must be a subset of: {', '.join(image_channels)}"

    def _run_patch(self, patch: Polygon) -> gpd.GeoDataFrame:
        """Run segmentation on one patch

        Args:
            patch: Patch, represented as a `shapely` polygon

        Returns:
            A list of cells, represented as polygons
        """
        bounds = [int(x) for x in patch.bounds]

        image = self.image.sel(
            c=self.channels,
            x=slice(bounds[0], bounds[2]),
            y=slice(bounds[1], bounds[3]),
        ).values

        assert np.issubdtype(
            image.dtype, np.integer
        ), f"Invalid image type {image.dtype}. Transform it to an integer dtype, e.g. `np.uint8`."

        if self.gaussian_sigma > 0:
            image = np.stack([gaussian_filter(c, sigma=self.gaussian_sigma) for c in image])
        if self.clip_limit > 0:
            image = np.stack(
                [
                    exposure.equalize_adapthist(
                        c,
                        clip_limit=self.clip_limit,
                        kernel_size=self.clahe_kernel_size,
                    )
                    for c in image
                ]
            )

        if patch.area < box(*bounds).area:
            image = image * shapes.rasterize(patch, image.shape[1:], bounds)

        cells = shapes.vectorize(self.method(image))
        cells.geometry = cells.translate(*bounds[:2])

        return cells[cells.area >= self.min_area] if self.min_area > 0 else cells

    def write_patch_cells(self, patch_dir: str, patch_index: int, recover: bool = False):
        """Run segmentation on one patch, and save the result in a dedicated directory

        Args:
            patch_dir: Directory inside which segmentation results will be saved
            patch_index: Index of the patch on which to run segmentation. NB: the number of patches is `len(sdata['image_patches'])`
            recover: If `True`, the function will not run segmentation if the output file already exists
        """
        patch_dir: Path = Path(patch_dir)
        patch_dir.mkdir(parents=True, exist_ok=True)
        output_path = patch_dir / f"{patch_index}.parquet"

        if recover and output_path.exists():
            return

        patch = self.patches_gdf.geometry[patch_index]

        cells = self._run_patch(patch)
        cells.to_parquet(output_path)

    def write_patches_cells(self, patch_dir: str, recover: bool = False):
        """Run segmentation on all patches, and save the result in a dedicated directory

        Args:
            patch_dir: Directory inside which segmentation results will be saved
            recover: If `True`, the function will not run segmentation on already-segmented patches
        """
        functions = [
            partial(self.write_patch_cells, patch_dir, patch_index, recover)
            for patch_index in range(len(self.patches_gdf))
        ]

        if settings.parallelization_backend is not None and not len(spatialdata.get_dask_backing_files(self.image)):
            log.warning(
                "Having backed data is recommended when running staining-based segmentation with a parallelization backend. "
                "Consider saving your data on disk with `sdata.write(...)`, and opening it with `sdata = spatialdata.read_zarr(...)`."
            )

        settings._run_with_backend(functions)

    @classmethod
    def read_patches_cells(cls, patch_dir: str | list[str]) -> gpd.GeoDataFrame:
        """Read all patch segmentation results after running `write_patch_cells` on all patches

        Args:
            patch_dir: Directory provided when running `write_patch_cells` containing the `.parquet` files. For multi-step segmentation, provide a list of directories (one per step).

        Returns:
            A list of cells represented as `shapely` polygons
        """
        patch_dirs = patch_dir if isinstance(patch_dir, list) else [patch_dir]
        files = [f for patch_dir in patch_dirs for f in Path(patch_dir).iterdir() if f.suffix == ".parquet"]

        cells = pd.concat([gpd.read_parquet(file) for file in files]).reset_index(drop=True)

        log.info(f"Found {len(cells)} total cells")
        return cells

    @classmethod
    def add_shapes(cls, sdata: SpatialData, cells: gpd.GeoDataFrame, image_key: str, key_added: str):
        """Adding `shapely` polygon to the `SpatialData` object

        Args:
            sdata: A `SpatialData` object
            cells: `GeoDataFrame` containing the cell boundaries after segmentation
            image_key: Key of the image on which segmentation has been run
            shapes_key: Name to provide to the geodataframe to be created
        """
        image = get_spatial_image(sdata, image_key)

        cells.index = image_key + cells.index.astype(str)

        cells = ShapesModel.parse(cells, transformations=get_transformation(image, get_all=True).copy())
        add_spatial_element(sdata, key_added, cells)

        log.info(f"Added {len(cells)} cell boundaries in sdata['{key_added}']")
