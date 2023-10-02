import numpy as np
import pandas as pd
from anndata import AnnData
from spatialdata import SpatialData
from spatialdata.models import TableModel

from .._constants import SopaKeys
from .._sdata import get_boundaries, get_intrinsic_cs, get_key, get_spatial_image
from . import shapes


def aggregate(
    sdata: SpatialData, gene_column: str | None, intensity_mean: bool = True, overwrite: bool = True
):
    image_key, image = get_spatial_image(sdata)
    shapes_key, geo_df = get_boundaries(sdata, return_key=True)

    table = table if sdata.table is None else sdata.table

    assert (
        intensity_mean or gene_column is not None or table is not None
    ), f"You must choose at least one aggregation: transcripts or fluorescence intensities"

    if gene_column is not None:
        points_key = get_key(sdata, "points")

        table = sdata.aggregate(
            values=points_key,
            by=shapes_key,
            value_key=gene_column,
            agg_func="count",
            target_coordinate_system=get_intrinsic_cs(sdata, points_key),
        ).table

    if intensity_mean:
        mean_intensities = shapes.average(image, geo_df.geometry)

    if table is None:
        table = AnnData(
            mean_intensities,
            dtype=mean_intensities.dtype,
            var=pd.DataFrame(index=image.c),
            obs=pd.DataFrame(index=geo_df.index),
        )
    elif intensity_mean:
        table.obsm[SopaKeys.INTENSITIES_OBSM] = pd.DataFrame(
            mean_intensities, columns=image.coords["c"].values, index=table.obs_names
        )

    table.obsm["spatial"] = np.array([[centroid.x, centroid.y] for centroid in geo_df.centroid])
    table.obs[SopaKeys.REGION_KEY] = pd.Series(shapes_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.SLIDE_KEY] = pd.Series(image_key, index=table.obs_names, dtype="category")
    table.obs[SopaKeys.INSTANCE_KEY] = geo_df.index

    if "spatialdata_attrs" in table.uns:
        del table.uns["spatialdata_attrs"]

    table = TableModel.parse(
        table,
        region_key=SopaKeys.REGION_KEY,
        region=shapes_key,
        instance_key=SopaKeys.INSTANCE_KEY,
    )

    if sdata.table is not None and overwrite:
        del sdata.table

    sdata.table = table
