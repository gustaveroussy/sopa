from __future__ import annotations

import logging
from pathlib import Path
from typing import Callable

import geopandas as gpd
import numpy as np
from scipy.ndimage import gaussian_filter
from shapely import affinity
from shapely.geometry import Polygon, box
from skimage import exposure
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from tqdm import tqdm

from .._constants import SopaKeys
from .._sdata import get_spatial_image, save_shapes
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
        gaussian_sigma: float = 1,
    ):
        """Generalized staining-based segmentation

        !!! note "Sequential segmentation (slower)"
            ```python
            from sopa.segmentation.stainings import StainingSegmentation

            method = ... # custom callable that runs segmentation on each patch

            segmentation = StainingSegmentation(sdata, method, "DAPI")
            segmentation.write_patches_cells("./temp_dir")

            cells = StainingSegmentation.read_patches_cells("./temp_dir")
            StainingSegmentation.add_shapes(sdata, cells, image_key, "method_name")
            ```

        !!! note "Parallel segmentation (faster)"
            ```python
            from sopa.segmentation.stainings import StainingSegmentation

            method = ... # custom callable that runs segmentation on each patch

            segmentation = StainingSegmentation(sdata, method, "DAPI")

            # Run all this in a parallel manner, e.g. on different jobs
            for i in range(len(sdata['sopa_patches'])):
                segmentation.write_patch_cells("./temp_dir", i)

            cells = StainingSegmentation.read_patches_cells("./temp_dir")
            StainingSegmentation.add_shapes(sdata, cells, image_key, "method_name")
            ```

        Args:
            sdata: A `SpatialData` object
            method: A segmentation `callable` whose input is an image of shape `(C, Y, X)` and output is a cell mask of shape `(Y, X)`. Each mask value `>0` represent a unique cell ID. Such callables can be found in `sopa.segmentation.methods`.
            channels: One or a list of channel names used for segmentation. If only one channel is provided, the image given to the `method` will be of shape `(1, Y, X)`.
            image_key: Optional key of `sdata` containing the image (no needed if there is only one image)
            min_area: Minimum area (in pixels^2) for a cell to be kept
            clip_limit: Parameter for skimage.exposure.equalize_adapthist (applied before running cellpose)
            gaussian_sigma: Parameter for scipy gaussian_filter (applied before running cellpose)
        """
        self.sdata = sdata
        self.method = method
        self.channels = [channels] if isinstance(channels, str) else channels

        self.min_area = min_area
        self.clip_limit = clip_limit
        self.gaussian_sigma = gaussian_sigma

        self.image_key, self.image = get_spatial_image(sdata, key=image_key, return_key=True)

        image_channels = self.image.coords["c"].values
        assert np.isin(
            channels, image_channels
        ).all(), f"Channel names must be a subset of: {', '.join(image_channels)}"

    def _run_patch(self, patch: Polygon) -> list[Polygon]:
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

        image = gaussian_filter(image, sigma=self.gaussian_sigma)
        image = exposure.equalize_adapthist(image, clip_limit=self.clip_limit)

        if patch.area < box(*bounds).area:
            image = image * shapes.rasterize(patch, image.shape[1:], bounds)

        cells = shapes.geometrize(self.method(image))

        if self.min_area > 0:
            cells = [cell for cell in cells if cell.area >= self.min_area]

        return [affinity.translate(cell, *bounds[:2]) for cell in cells]

    def write_patch_cells(self, patch_dir: str, patch_index: int):
        """Run segmentation on one patch, and save the result in a dedicated directory

        Args:
            patch_dir: Directory inside which segmentation results will be saved
            patch_index: Index of the patch on which to run segmentation. NB: the number of patches is `len(sdata['sopa_patches'])`
        """
        patch = self.sdata[SopaKeys.PATCHES].geometry[patch_index]

        cells = self._run_patch(patch)
        gdf = gpd.GeoDataFrame(geometry=cells)

        patch_dir: Path = Path(patch_dir)
        patch_dir.mkdir(parents=True, exist_ok=True)
        patch_file = patch_dir / f"{patch_index}.parquet"

        gdf.to_parquet(patch_file)

    def write_patches_cells(self, patch_dir: str):
        log.warn(
            "Running segmentation in a sequential manner. This is not recommended on large images because it can be extremely slow (see https://github.com/gustaveroussy/sopa/discussions/36 for more details)"
        )
        for patch_index in tqdm(range(len(self.sdata[SopaKeys.PATCHES])), desc="Run all patches"):
            self.write_patch_cells(patch_dir, patch_index)

    @classmethod
    def read_patches_cells(cls, patch_dir: str | list[str]) -> list[Polygon]:
        """Read all patch segmentation results after running `write_patch_cells` on all patches

        Args:
            patch_dir: Directory provided when running `write_patch_cells` containing the `.parquet` files. For multi-step segmentation, provide a list of directories (one per step).

        Returns:
            A list of cells represented as `shapely` polygons
        """
        cells = []

        files = [f for f in Path(patch_dir).iterdir() if f.suffix == ".parquet"]
        for file in tqdm(files, desc="Reading patches"):
            cells += list(gpd.read_parquet(file).geometry)

        log.info(f"Found {len(cells)} total cells")
        return cells

    @classmethod
    def add_shapes(cls, sdata: SpatialData, cells: list[Polygon], image_key: str, shapes_key: str):
        """Adding `shapely` polygon to the `SpatialData` object

        Args:
            sdata: A `SpatialData` object
            cells: List of polygons after segmentation
            image_key: Key of the image on which segmentation has been run
            shapes_key: Name to provide to the geodataframe to be created
        """
        image = get_spatial_image(sdata, image_key)

        geo_df = gpd.GeoDataFrame({"geometry": cells})
        geo_df.index = image_key + geo_df.index.astype(str)

        geo_df = ShapesModel.parse(
            geo_df, transformations=get_transformation(image, get_all=True).copy()
        )
        sdata.shapes[shapes_key] = geo_df
        save_shapes(sdata, shapes_key, overwrite=True)

        log.info(f"Added {len(geo_df)} cell boundaries in sdata['{shapes_key}']")
