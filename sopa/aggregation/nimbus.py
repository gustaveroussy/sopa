import logging
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from nimbus_inference.nimbus import Nimbus, prep_naming_convention
from nimbus_inference.utils import MultiplexDataset
from shapely.geometry import Polygon
from spatialdata import SpatialData
from tifffile import imwrite
from xarray import DataArray

from ..constants import NimbusFiles
from ..shapes import expand_radius, rasterize_labeled
from ..utils import get_boundaries, get_cache_dir, get_spatial_image, to_intrinsic, validated_channel_names

log = logging.getLogger(__name__)


def nimbus_aggregation(
    sdata: SpatialData,
    image_key: str | None = None,
    shapes_key: str | None = None,
    include_channels: list[str] | None = None,
    expand_radius_ratio: float = 0,
    no_overlap: bool = False,
    quantile: float = 0.999,
    multiprocessing: bool = True,
    batch_size: int = 4,
    clip_values: tuple = (0, 2),
    delete_cache: bool = True,
    **nimbus_model_kwargs,
) -> AnnData:
    """Run [Nimbus](https://nimbus-inference.readthedocs.io/en/latest/?badge=latest) aggregation on a SpatialData object, and predict marker confidence scores for each cell.

    !!! warning "Nimbus installation"
        Make sure to install the Nimbus-Inference in your environment (`pip install 'Nimbus-Inference'`) for this method to work.

    Args:
        sdata: A `SpatialData` object
        image_key: Key of `sdata` containing the image. If only one `images` element, this does not have to be provided.
        shapes_key: Key of `sdata` containing the cell boundaries. If only one `shapes` element, this does not have to be provided.
        include_channels: List of channels to include in the prediction. If None, all channels are included.
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate boundary stainings.
        no_overlap: If `True`, the (expanded) cells will not overlap.
        quantile: Quantile used for Nimbus intensity normalization.
        multiprocessing: Whether to use multiprocessing.
        batch_size: Nimbus inference batch size.
        clip_values: Values to clip images to after Nimbus intensity normalization.
        delete_cache: Whether to delete temporary cache files after inference.


    Returns:
        An `AnnData` object with the predicted marker confidence scores
    """
    image = get_spatial_image(sdata, image_key)

    geo_df = get_boundaries(sdata, key=shapes_key)
    geo_df = to_intrinsic(sdata, geo_df, image)
    geo_df = expand_radius(geo_df, expand_radius_ratio, no_overlap=no_overlap)
    cache_dir = get_cache_dir(sdata)

    return _run_nimbus(
        image=image,
        geo_df=geo_df,
        include_channels=include_channels,
        work_dir=cache_dir,
        quantile=quantile,
        batch_size=batch_size,
        clip_values=clip_values,
        multiprocessing=multiprocessing,
        delete_cache=delete_cache,
        **nimbus_model_kwargs,
    )


def _write_segmentation_data(image: DataArray, mask: np.ndarray, work_dir: Path, channel_names: list[str]) -> None:
    tiff_dir = work_dir / NimbusFiles.IMAGE_DIR
    seg_dir = work_dir / NimbusFiles.SEGMENTATION_DIR
    tiff_dir.mkdir(parents=True, exist_ok=True)
    seg_dir.mkdir(parents=True, exist_ok=True)
    for stale_tiff in tiff_dir.glob("*.tiff"):
        stale_tiff.unlink()

    for c, name in enumerate(channel_names):
        imwrite(tiff_dir / f"{name}.tiff", image[c])

    seg_file = seg_dir / NimbusFiles.SEGMENTATION_MASK
    if seg_file.exists():
        seg_file.unlink()
    imwrite(seg_file, mask)


def _build_multiplex_dataset(work_dir: Path, include_channels: list[str] | None) -> MultiplexDataset:
    tiff_dir = work_dir / "image_data"
    seg_dir = work_dir / NimbusFiles.SEGMENTATION_DIR
    nimbus_output_dir = work_dir / "nimbus_output"
    nimbus_output_dir.mkdir(parents=True, exist_ok=True)

    fov_paths = sorted(str(p) for p in tiff_dir.iterdir() if p.is_dir() and not p.name.startswith("."))
    if not fov_paths:
        raise FileNotFoundError(f"No FOV directories found in {tiff_dir}")

    segmentation_naming = prep_naming_convention(seg_dir)

    missing = [p for p in fov_paths if not Path(segmentation_naming(str(p))).exists()]
    if missing:
        raise FileNotFoundError(f"Missing segmentation files for {len(missing)} FOV(s), first: {missing[0]}")

    return MultiplexDataset(
        fov_paths=[str(p) for p in fov_paths],
        suffix=".tiff",
        include_channels=include_channels,
        segmentation_naming_convention=segmentation_naming,
        output_dir=nimbus_output_dir,
    )


def _nimbus_inference(
    work_dir: Path,
    mpx_data: MultiplexDataset,
    quantile: float,
    batch_size: int,
    clip_values: tuple,
    multiprocessing: bool,
    **nimbus_model_kwargs,
) -> pd.DataFrame:
    nimbus_output_dir = work_dir / "nimbus_output"

    nimbus = Nimbus(
        dataset=mpx_data,
        save_predictions=False,
        batch_size=batch_size,
        output_dir=nimbus_output_dir,
        **nimbus_model_kwargs,
    )

    mpx_data.prepare_normalization_dict(
        quantile=quantile, clip_values=clip_values, multiprocessing=multiprocessing, n_subset=1, overwrite=True
    )

    return nimbus.predict_fovs()


def _run_nimbus(
    image: DataArray,
    geo_df: gpd.GeoDataFrame | list[Polygon],
    include_channels: list[str] | None,
    work_dir: Path,
    quantile: float,
    batch_size: int,
    clip_values: tuple,
    multiprocessing: bool,
    delete_cache: bool = True,
    **nimbus_model_kwargs,
) -> AnnData:
    available_channels = validated_channel_names(image)
    channel_names = available_channels if include_channels is None else include_channels
    unknown_channels = sorted(set(channel_names) - set(available_channels))
    if unknown_channels:
        raise ValueError(f"Unknown include_channels: {unknown_channels}. Available channels: {available_channels}")

    log.info(f"Predicting marker confidence scores for {len(geo_df)} cells with {include_channels=}")

    if image.shape[1] < 2 or image.shape[2] < 2:
        raise ValueError(f"Image spatial dimensions must be >= 2, got {(image.shape[1], image.shape[2])}.")

    labels = np.arange(1, len(geo_df) + 1, dtype=np.int32)

    mask = rasterize_labeled(
        shapes=((geom, label) for geom, label in zip(geo_df.geometry, labels)),
        out_shape=image.shape[1:],
    )

    _write_segmentation_data(image, mask, work_dir, channel_names)

    mpx_data = _build_multiplex_dataset(work_dir, channel_names)

    mpx_data.check_inputs()

    cell_table = _nimbus_inference(
        work_dir,
        mpx_data,
        quantile=quantile,
        batch_size=batch_size,
        clip_values=clip_values,
        multiprocessing=multiprocessing,
        **nimbus_model_kwargs,
    )

    if delete_cache and work_dir.exists():
        shutil.rmtree(work_dir)

    X = cell_table[channel_names].to_numpy()

    adata = AnnData(
        X,
        var=pd.DataFrame(index=channel_names),
        obs=pd.DataFrame(index=geo_df.index.astype(str)),
    )

    return adata
