import logging
from pathlib import Path
from typing import Callable

import geopandas as gpd
import numpy as np
import zarr
from scipy.ndimage import gaussian_filter
from shapely import affinity
from shapely.geometry import Polygon, box
from skimage import exposure
from spatialdata import SpatialData
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from tqdm import tqdm

from .._constants import SopaKeys
from .._sdata import get_spatial_image
from . import shapes
from .patching import Patches2D

log = logging.getLogger(__name__)


class StainingSegmentation:
    def __init__(
        self,
        sdata: SpatialData,
        method: Callable,
        channels: list[str] | str,
        min_area: float = 0,
    ):
        """Generalized staining-based segmentation

        !!! note "Sequential segmentation (slower)"
            ```python
            from sopa.segmentation.stainings import StainingSegmentation

            method = ... # custom callable that runs segmentation on each patch

            segmentation = StainingSegmentation(sdata, method, "DAPI")
            cells = segmentation.run_patches(2000, 100)
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
            channels: One or a list of channel names used for segmentation
            min_area: Minimum area (in pixels^2) for a cell to be kept
        """
        self.sdata = sdata
        self.method = method
        self.channels = [channels] if isinstance(channels, str) else channels
        self.min_area = min_area

        self.image_key, self.image = get_spatial_image(sdata, return_key=True)

        image_channels = self.image.coords["c"].values
        assert np.isin(
            channels, image_channels
        ).all(), f"Channel names must be a subset of: {', '.join(image_channels)}"

    def _run_patch(
        self, patch: Polygon, clip_limit: float = 0.2, sigma: float = 1
    ) -> list[Polygon]:
        """Run segmentation on one patch

        Args:
            patch: Patch, represented as a `shapely` polygon
            clip_limit: parameter for skimage.exposure.equalize_adapthist
            sigma: parameter for scipy gaussian_filter

        Returns:
            A list of cells, represented as polygons
        """
        bounds = [int(x) for x in patch.bounds]

        image = self.image.sel(
            c=self.channels,
            x=slice(bounds[0], bounds[2]),
            y=slice(bounds[1], bounds[3]),
        ).values

        image = gaussian_filter(image, sigma=sigma)
        image = exposure.equalize_adapthist(image, clip_limit=clip_limit)

        if patch.area < box(*bounds).area:
            image = image * shapes.rasterize(patch, image.shape[1:], bounds)

        cells = shapes.geometrize(self.method(image))
        cells = shapes.filter(cells, self.min_area)

        return [affinity.translate(cell, *bounds[:2]) for cell in cells]

    def write_patch_cells(self, patch_dir: str, patch_index: int):
        """Run segmentation on one patch, and save the result in a dedicated directory

        Args:
            patch_dir: Directory inside which segmentation results will be saved
            patch_index: Index of the patch on which to run segmentation. NB: the number of patches is `len(sdata['sopa_patches'])`
        """
        patch = self.sdata[SopaKeys.PATCHES].geometry[patch_index]
        cells = self._run_patch(patch)

        patch_file = Path(patch_dir) / f"{patch_index}.zarr.zip"

        with zarr.ZipStore(patch_file, mode="w") as store:
            g = zarr.group(store=store)

            for i, cell in enumerate(cells):
                coords = np.array(cell.exterior.coords)
                g.array(f"cell_{i}", coords, dtype=coords.dtype, chunks=coords.shape)

    def run_patches(
        self,
        patch_width: int,
        patch_overlap: int,
    ) -> list[Polygon]:
        """Run segmentation over all patches, in a sequential manner (this is slower than running all patches in parallel)

        Args:
            patch_width: Width of the patches
            patch_overlap: Number of pixels of overlap between patches

        Returns:
            A list of cells represented as `shapely` polygons
        """
        self.patches = Patches2D(self.sdata, self.image_key, patch_width, patch_overlap)

        cells = [
            cell
            for patch in tqdm(self.patches.polygons, desc="Run on patches")
            for cell in self._run_patch(patch)
        ]
        cells = shapes.solve_conflicts(cells)
        return cells

    @classmethod
    def read_patches_cells(cls, patch_dir: str) -> list[Polygon]:
        """Read all patch segmentation results after running `write_patch_cells` on all patches

        Args:
            patch_dir: Directory provided when running `write_patch_cells`

        Returns:
            A list of cells represented as `shapely` polygons
        """
        cells = []

        files = [f for f in Path(patch_dir).iterdir() if f.suffix == ".zip"]
        for file in tqdm(files, desc="Reading patches"):
            z = zarr.open(file, mode="r")
            for _, coords_zarr in z.arrays():
                cells.append(Polygon(coords_zarr[:]))

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

        geo_df = ShapesModel.parse(geo_df, transformations=get_transformation(image, get_all=True))
        sdata.add_shapes(shapes_key, geo_df, overwrite=True)

        log.info(f"Added {len(geo_df)} cell boundaries in sdata['{shapes_key}']")
