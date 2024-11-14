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
):
    """Overlay a segmentation on top of an existing segmentation

    Args:
        sdata: A `SpatialData` object
        shapes_key: The key of the new shapes to be added
        gene_column: Key of the points dataframe containing the genes names
        area_ratio_threshold: Threshold between 0 and 1. For each original cell overlapping with a new cell, we compute the overlap-area/cell-area, if above the threshold the cell is removed.
        image_key: Optional key of the original image
    """
    aggregate_genes, aggregate_channels = False, False

    assert SopaKeys.TABLE in sdata.tables, "No table found in the SpatialData object"

    assert SopaKeys.UNS_KEY in sdata.tables["table"].uns, "It seems the table was not aggregated using `sopa.aggregate`"

    sopa_attrs = sdata.tables["table"].uns[SopaKeys.UNS_KEY]

    aggregate_genes = sopa_attrs[SopaKeys.UNS_HAS_TRANSCRIPTS]
    aggregate_channels = sopa_attrs[SopaKeys.UNS_HAS_INTENSITIES]

    if aggregate_genes and gene_column is None:
        points = get_spatial_element(sdata.points, key=sdata.attrs.get(SopaAttrs.TRANSCRIPTS))
        gene_column = get_feature_key(points, raise_error=True)

    aggr = Aggregator(sdata, image_key=image_key, shapes_key=shapes_key)

    old_table: AnnData = aggr.sdata.tables[SopaKeys.TABLE]
    aggr.sdata.tables[SopaKeys.OLD_TABLE] = old_table
    del aggr.sdata.tables[SopaKeys.TABLE]

    old_shapes_key = old_table.uns["spatialdata_attrs"]["region"]
    instance_key = old_table.uns["spatialdata_attrs"]["instance_key"]

    if isinstance(old_shapes_key, list):
        assert len(old_shapes_key) == 1, "Can't overlap segmentation on multi-region SpatialData object"
        old_shapes_key = old_shapes_key[0]

    old_geo_df = aggr.sdata[old_shapes_key]
    geo_df = to_intrinsic(aggr.sdata, aggr.geo_df, old_geo_df)

    geo_df.index.name = None
    gdf_join = gpd.sjoin(old_geo_df, geo_df)
    gdf_join["geometry_right"] = gdf_join["index_right"].map(lambda i: geo_df.geometry.iloc[i])
    gdf_join["overlap_ratio"] = gdf_join.apply(_overlap_area_ratio, axis=1)
    gdf_join: gpd.GeoDataFrame = gdf_join[gdf_join.overlap_ratio >= area_ratio_threshold]

    table_crop = old_table[~np.isin(old_table.obs[instance_key], gdf_join.index)].copy()
    table_crop.obs[SopaKeys.CELL_OVERLAY_KEY] = False

    aggr.compute_table(aggregate_channels=aggregate_channels, aggregate_genes=aggregate_genes, gene_column=gene_column)
    aggr.table.obs[SopaKeys.CELL_OVERLAY_KEY] = True

    aggr.table = anndata.concat([table_crop, aggr.table], uns_merge="first", join="outer", fill_value=0)
    _fillna(aggr.table.obs)

    aggr.shapes_key = f"{old_shapes_key}_overlay_{aggr.shapes_key}"
    geo_df_cropped = old_geo_df.loc[~old_geo_df.index.isin(gdf_join.index)]
    aggr.geo_df = pd.concat([geo_df_cropped, geo_df], join="outer", axis=0)
    aggr.geo_df.attrs = old_geo_df.attrs

    aggr.add_standardized_table()


def _overlap_area_ratio(row) -> float:
    poly: Polygon = row["geometry"]
    poly_right: Polygon = row["geometry_right"]
    return poly.intersection(poly_right).area / poly.area


def _fillna(df: pd.DataFrame):
    for key in df:
        if df[key].dtype == "category":
            df[key] = df[key].cat.add_categories("NA").fillna("NA")
        else:
            df[key] = df[key].fillna(0)
