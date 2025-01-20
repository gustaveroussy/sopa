import logging

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from shapely import Polygon
from spatialdata import SpatialData

from .._constants import SopaAttrs, SopaKeys
from ..utils import get_feature_key, get_spatial_element, to_intrinsic
from . import Aggregator

log = logging.getLogger(__name__)


def overlay_segmentation(
    sdata: SpatialData,
    shapes_key: str,
    gene_column: str | None = None,
    area_ratio_threshold: float = 0.25,
    image_key: str | None = None,
    table_key: str = SopaKeys.TABLE,
):
    """Overlay a segmentation on top of an existing segmentation

    Args:
        sdata: A `SpatialData` object
        shapes_key: The key of the new shapes to be added
        gene_column: Key of the points dataframe containing the genes names
        area_ratio_threshold: Threshold between 0 and 1. For each original cell overlapping with a new cell, we compute the overlap-area/cell-area, if above the threshold the cell is removed.
        image_key: Optional key of the original image
        table_key: Key of the table to be overlayed
    """
    aggregate_genes, aggregate_channels = False, False

    assert table_key in sdata.tables, f"No table with name '{table_key}' found in the SpatialData object"

    old_table: AnnData = sdata.tables[table_key]

    assert SopaKeys.UNS_KEY in old_table.uns, "It seems the table was not aggregated using `sopa.aggregate`"

    sopa_attrs = old_table.uns[SopaKeys.UNS_KEY]

    aggregate_genes = sopa_attrs[SopaKeys.UNS_HAS_TRANSCRIPTS]
    aggregate_channels = sopa_attrs[SopaKeys.UNS_HAS_INTENSITIES]

    if aggregate_genes and gene_column is None:
        points = get_spatial_element(sdata.points, key=sdata.attrs.get(SopaAttrs.TRANSCRIPTS))
        gene_column = get_feature_key(points, raise_error=True)

    aggregator = Aggregator(sdata, image_key=image_key, shapes_key=shapes_key)
    aggregator.sdata.tables[f"{SopaKeys.OLD_TABLE_PREFFIX}{table_key}"] = old_table
    del aggregator.sdata.tables[table_key]

    old_shapes_key = old_table.uns["spatialdata_attrs"]["region"]
    instance_key = old_table.uns["spatialdata_attrs"]["instance_key"]

    if isinstance(old_shapes_key, list):
        assert len(old_shapes_key) == 1, "Can't overlap segmentation on multi-region SpatialData object"
        old_shapes_key = old_shapes_key[0]

    old_geo_df = aggregator.sdata[old_shapes_key]
    geo_df = to_intrinsic(aggregator.sdata, aggregator.geo_df, old_geo_df)

    geo_df.index.name = None
    gdf_join = gpd.sjoin(old_geo_df, geo_df)
    gdf_join["geometry_right"] = gdf_join["index_right"].map(lambda i: geo_df.geometry.iloc[i])
    gdf_join["overlap_ratio"] = gdf_join.apply(_overlap_area_ratio, axis=1)
    gdf_join: gpd.GeoDataFrame = gdf_join[gdf_join.overlap_ratio >= area_ratio_threshold]

    table_crop = old_table[~np.isin(old_table.obs[instance_key], gdf_join.index)].copy()
    table_crop.obs[SopaKeys.CELL_OVERLAY_KEY] = False

    aggregator.compute_table(
        aggregate_channels=aggregate_channels,
        aggregate_genes=aggregate_genes,
        gene_column=gene_column,
        key_added=table_key,
    )
    aggregator.table.obs[SopaKeys.CELL_OVERLAY_KEY] = True

    aggregator.table = anndata.concat(
        [table_crop, aggregator.table],
        uns_merge="first",
        join="outer",
    )

    aggregator.shapes_key = f"{old_shapes_key}_overlay_{aggregator.shapes_key}"

    geo_df_cropped = old_geo_df.loc[~old_geo_df.index.isin(gdf_join.index)]
    aggregator.geo_df = pd.concat([geo_df_cropped, geo_df], join="outer", axis=0)
    aggregator.geo_df.attrs = old_geo_df.attrs

    aggregator.add_standardized_table(table_key)


def _overlap_area_ratio(row) -> float:
    poly: Polygon = row["geometry"]
    poly_right: Polygon = row["geometry_right"]
    return poly.intersection(poly_right).area / poly.area
