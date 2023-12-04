# Readers for spatial-transcriptomics technologies
# Updated from spatialdata-io: https://spatialdata.scverse.org/projects/io/en/latest/
# In the future, we will completely rely on spatialdata-io (when stable enough)

import json
import re
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import Any

import anndata
import dask.dataframe as dd
import geopandas
import numpy as np
import pandas as pd
import spatialdata_io
from dask import array as da
from dask.dataframe import read_parquet
from dask_image.imread import imread
from spatialdata import SpatialData
from spatialdata._logging import logger
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel
from spatialdata.transformations import Affine, Identity
from spatialdata.transformations.transformations import Identity, Scale
from spatialdata_io._constants._constants import MerscopeKeys, XeniumKeys


def _get_channel_names(images_dir: Path) -> list[str]:
    exp = r"mosaic_(?P<stain>[\w|-]+[0-9]?)_z(?P<z>[0-9]+).tif"
    matches = [re.search(exp, file.name) for file in images_dir.iterdir()]

    stainings = {match.group("stain") for match in matches if match}

    return list(stainings)


def _get_file_paths(
    path: Path, vpt_outputs: Path | str | dict[str, Any] | None
) -> tuple[Path, Path, Path]:
    """
    Gets the MERSCOPE file paths when vpt_outputs is provided

    That is, (i) the file of transcript per cell, (ii) the cell metadata file, and (iii) the cell boundary file
    """
    if vpt_outputs is None:
        return (
            path / MerscopeKeys.COUNTS_FILE,
            path / MerscopeKeys.CELL_METADATA_FILE,
            path / MerscopeKeys.BOUNDARIES_FILE,
        )

    if isinstance(vpt_outputs, str) or isinstance(vpt_outputs, Path):
        vpt_outputs = Path(vpt_outputs)

        plausible_boundaries = [
            vpt_outputs / MerscopeKeys.CELLPOSE_BOUNDARIES,
            vpt_outputs / MerscopeKeys.WATERSHED_BOUNDARIES,
        ]
        valid_boundaries = [path for path in plausible_boundaries if path.exists()]

        assert (
            valid_boundaries
        ), f"Boundary file not found - expected to find one of these files: {', '.join(map(str, plausible_boundaries))}"

        return (
            vpt_outputs / MerscopeKeys.COUNTS_FILE,
            vpt_outputs / MerscopeKeys.CELL_METADATA_FILE,
            valid_boundaries[0],
        )

    if isinstance(vpt_outputs, dict):
        return (
            vpt_outputs[MerscopeKeys.VPT_NAME_COUNTS],
            vpt_outputs[MerscopeKeys.VPT_NAME_OBS],
            vpt_outputs[MerscopeKeys.VPT_NAME_BOUNDARIES],
        )

    raise ValueError(
        f"`vpt_outputs` has to be either `None`, `str`, `Path`, or `dict`. Found type {type(vpt_outputs)}."
    )


def merscope(
    path: str | Path,
    vpt_outputs: Path | str | dict[str, Any] | None = None,
    z_layers: int | list[int] | None = 3,
    region_name: str | None = None,
    slide_name: str | None = None,
    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> SpatialData:
    """Read MERSCOPE data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.merscope.html).

    Args:
        path: Path to the MERSCOPE directory containing all the experiment files
        **kwargs: See link above.

    Returns:
        A `SpatialData` object representing the MERSCOPE experiment
    """
    if "chunks" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["chunks"] = (1, 4096, 4096)
    if "scale_factors" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    path = Path(path).absolute()
    count_path, obs_path, boundaries_path = _get_file_paths(path, vpt_outputs)
    images_dir = path / MerscopeKeys.IMAGES_DIR

    microns_to_pixels = Affine(
        np.genfromtxt(images_dir / MerscopeKeys.TRANSFORMATION_FILE),
        input_axes=("x", "y"),
        output_axes=("x", "y"),
    )

    vizgen_region = path.name if region_name is None else region_name
    slide_name = path.parent.name if slide_name is None else slide_name
    dataset_id = f"{slide_name}_{vizgen_region}"
    region = f"{dataset_id}_polygons"

    # Images
    images = {}

    z_layers = [z_layers] if isinstance(z_layers, int) else z_layers or []

    stainings = _get_channel_names(images_dir)
    if stainings:
        for z_layer in z_layers:
            im = da.stack(
                [
                    imread(images_dir / f"mosaic_{stain}_z{z_layer}.tif", **imread_kwargs).squeeze()
                    for stain in stainings
                ],
                axis=0,
            )
            parsed_im = Image2DModel.parse(
                im,
                dims=("c", "y", "x"),
                transformations={"microns": microns_to_pixels.inverse()},
                c_coords=stainings,
                **image_models_kwargs,
            )
            images[f"{dataset_id}_z{z_layer}"] = parsed_im

    # Transcripts
    points = {}
    transcript_path = path / MerscopeKeys.TRANSCRIPTS_FILE
    if transcript_path.exists():
        points[f"{dataset_id}_transcripts"] = _get_points(transcript_path)
    else:
        logger.warning(
            f"Transcript file {transcript_path} does not exist. Transcripts are not loaded."
        )

    # Polygons
    shapes = {}
    if boundaries_path.exists():
        shapes[f"{dataset_id}_polygons"] = _get_polygons(boundaries_path)
    else:
        logger.warning(
            f"Boundary file {boundaries_path} does not exist. Cell boundaries are not loaded."
        )

    # Table
    table = None
    if count_path.exists() and obs_path.exists():
        table = _get_table(count_path, obs_path, vizgen_region, slide_name, dataset_id, region)
    else:
        logger.warning(
            f"At least one of the following files does not exist: {count_path}, {obs_path}. The table is not loaded."
        )

    return SpatialData(shapes=shapes, points=points, images=images, table=table)


def _get_points(transcript_path: Path):
    transcript_df = dd.read_csv(transcript_path)
    transcripts = PointsModel.parse(
        transcript_df,
        coordinates={"x": MerscopeKeys.GLOBAL_X, "y": MerscopeKeys.GLOBAL_Y},
        transformations={"microns": Identity()},
    )
    transcripts["gene"] = transcripts["gene"].astype("category")
    return transcripts


def _get_polygons(boundaries_path: Path) -> geopandas.GeoDataFrame:
    geo_df = geopandas.read_parquet(boundaries_path)
    geo_df = geo_df.rename_geometry("geometry")
    geo_df = geo_df[geo_df[MerscopeKeys.Z_INDEX] == 0]  # Avoid duplicate boundaries on all z-levels
    geo_df.geometry = geo_df.geometry.map(
        lambda x: x.geoms[0]
    )  # The MultiPolygons contain only one polygon
    geo_df.index = geo_df[MerscopeKeys.METADATA_CELL_KEY].astype(str)

    return ShapesModel.parse(geo_df, transformations={"microns": Identity()})


def _get_table(
    count_path: Path,
    obs_path: Path,
    vizgen_region: str,
    slide_name: str,
    dataset_id: str,
    region: str,
) -> anndata.AnnData:
    data = pd.read_csv(count_path, index_col=0, dtype={MerscopeKeys.COUNTS_CELL_KEY: str})
    obs = pd.read_csv(obs_path, index_col=0, dtype={MerscopeKeys.METADATA_CELL_KEY: str})

    is_gene = ~data.columns.str.lower().str.contains("blank")
    adata = anndata.AnnData(data.loc[:, is_gene], dtype=data.values.dtype, obs=obs)

    adata.obsm["blank"] = data.loc[:, ~is_gene]  # blank fields are excluded from adata.X
    adata.obsm["spatial"] = adata.obs[[MerscopeKeys.CELL_X, MerscopeKeys.CELL_Y]].values
    adata.obs["region"] = pd.Series(vizgen_region, index=adata.obs_names, dtype="category")
    adata.obs["slide"] = pd.Series(slide_name, index=adata.obs_names, dtype="category")
    adata.obs["dataset_id"] = pd.Series(dataset_id, index=adata.obs_names, dtype="category")
    adata.obs[MerscopeKeys.REGION_KEY] = pd.Series(region, index=adata.obs_names, dtype="category")
    adata.obs[MerscopeKeys.METADATA_CELL_KEY] = adata.obs.index

    table = TableModel.parse(
        adata,
        region_key=MerscopeKeys.REGION_KEY.value,
        region=region,
        instance_key=MerscopeKeys.METADATA_CELL_KEY.value,
    )
    return table


def xenium(
    path: str | Path,
    imread_kwargs=MappingProxyType({}),
    image_models_kwargs=MappingProxyType({}),
) -> SpatialData:
    """Read Xenium data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.xenium.html).

    Args:
        path: Path to the Xenium directory containing all the experiment files
        imread_kwargs: See link above.
        image_models_kwargs:See link above.

    Returns:
        A `SpatialData` object representing the Xenium experiment
    """
    if "chunks" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["chunks"] = (1, 4096, 4096)
    if "scale_factors" not in image_models_kwargs:
        if isinstance(image_models_kwargs, MappingProxyType):
            image_models_kwargs = {}
        assert isinstance(image_models_kwargs, dict)
        image_models_kwargs["scale_factors"] = [2, 2, 2, 2]

    path = Path(path)
    with open(path / XeniumKeys.XENIUM_SPECS) as f:
        specs = json.load(f)

    points = {"transcripts": _get_points_xenium(path, specs)}

    images = {
        "morphology_mip": _get_images_xenium(
            path,
            XeniumKeys.MORPHOLOGY_MIP_FILE,
            imread_kwargs,
            image_models_kwargs,
        )
    }

    return SpatialData(images=images, points=points)


def _get_points_xenium(path: Path, specs: dict[str, Any]):
    table = read_parquet(path / XeniumKeys.TRANSCRIPTS_FILE)
    table["feature_name"] = table["feature_name"].apply(
        lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x),
        meta=("feature_name", "object"),
    )

    transform = Scale([1.0 / specs["pixel_size"], 1.0 / specs["pixel_size"]], axes=("x", "y"))
    points = PointsModel.parse(
        table,
        coordinates={
            "x": XeniumKeys.TRANSCRIPTS_X,
            "y": XeniumKeys.TRANSCRIPTS_Y,
            "z": XeniumKeys.TRANSCRIPTS_Z,
        },
        feature_key=XeniumKeys.FEATURE_NAME,
        instance_key=XeniumKeys.CELL_ID,
        transformations={"global": transform},
    )
    return points


def _get_images_xenium(
    path: Path,
    file: str,
    imread_kwargs: Mapping[str, Any] = MappingProxyType({}),
    image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    image = imread(path / file, **imread_kwargs)
    return Image2DModel.parse(
        image,
        transformations={"global": Identity()},
        dims=("c", "y", "x"),
        c_coords=list(map(str, range(len(image)))),
        **image_models_kwargs,
    )


def cosmx(path: str, **kwargs: int) -> SpatialData:
    """Alias to the [spatialdata-io reader](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.cosmx.html).

    Args:
        path: Path to the CosMX data directory
        **kwargs: See link above.

    Returns:
        A `SpatialData` object representing the CosMX experiment
    """
    return spatialdata_io.cosmx(path, **kwargs)
