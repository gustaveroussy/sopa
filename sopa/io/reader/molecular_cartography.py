from pathlib import Path

import dask.array as da
import pandas as pd
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel

from ...constants import SopaAttrs
from .utils import _default_image_kwargs


def molecular_cartography(
    path: str | Path,
    region: str,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
) -> SpatialData:
    """Read *Molecular Cartography* data from *Resolve Bioscience* as a `SpatialData` object.

    Args:
        path: Path to the directory containing the `.tiff` images and `_results.txt` files.
        region: Name of the region to read. The region name can be found before the `_results.txt` file, e.g. `A2-1`.
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.

    Returns:
        A `SpatialData` object representing the *Resolve Bioscience* experiment
    """
    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    path = Path(path)
    dataset_id = _get_dataset_id(path, region)

    # Read the points
    transcripts = pd.read_csv(path / f"{dataset_id}_results.txt", sep="\t", header=None)
    transcripts.columns = ["x", "y", "z", "target_name", "unnamed"]

    transcripts = PointsModel.parse(transcripts, feature_key="target_name")
    transcripts_name = f"{dataset_id}_points"

    # Read the images
    images_paths = list(path.glob(f"{dataset_id}_*.tiff"))
    c_coords = [image_path.stem.split("_")[-1] for image_path in images_paths]

    image = Image2DModel.parse(
        da.concatenate([imread(image_path, **imread_kwargs) for image_path in images_paths], axis=0),
        dims=("c", "y", "x"),
        c_coords=c_coords,
        rgb=None,
        **image_models_kwargs,
    )
    image_name = f"{dataset_id}_image"

    return SpatialData(
        images={image_name: image},
        points={transcripts_name: transcripts},
        attrs={
            SopaAttrs.CELL_SEGMENTATION: image_name,
            SopaAttrs.TRANSCRIPTS: transcripts_name,
        },
    )


def _get_dataset_id(path: Path, region: str) -> str:
    _dataset_ids = [path.name[:-12] for path in path.glob("*_results.txt")]
    region_to_id = {dataset_id.split("_")[-1]: dataset_id for dataset_id in _dataset_ids}

    assert region in region_to_id, f"Region {region} not found. Must be one of {list(region_to_id.keys())}"

    return region_to_id[region]
