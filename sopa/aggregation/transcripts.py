import logging
from functools import partial

import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
from anndata import AnnData
from dask.diagnostics import ProgressBar
from scipy.sparse import coo_matrix
from spatialdata import SpatialData

from .. import settings
from .._constants import SopaAttrs
from ..utils import get_boundaries, get_feature_key, get_spatial_element, to_intrinsic

log = logging.getLogger(__name__)


def count_transcripts(
    sdata: SpatialData,
    gene_column: str | None = None,
    shapes_key: str | None = None,
    points_key: str | None = None,
    geo_df: gpd.GeoDataFrame | None = None,
) -> AnnData:
    """Counts transcripts per cell.

    Args:
        sdata: A `SpatialData` object
        gene_column: Column of the transcript dataframe containing the gene names
        shapes_key: Key of `sdata` containing the cell boundaries. If only one `shapes` element, this does not have to be provided.
        points_key: Key of `sdata` containing the transcripts. If only one `points` element, this does not have to be provided.
        geo_df: If the cell boundaries are not yet in `sdata`, a `GeoDataFrame` can be directly provided for cell boundaries

    Returns:
        An `AnnData` object of shape `(n_cells, n_genes)` with the counts per cell
    """
    points_key, points = get_spatial_element(
        sdata.points, key=points_key or sdata.attrs.get(SopaAttrs.TRANSCRIPTS), return_key=True
    )

    if geo_df is None:
        geo_df = get_boundaries(sdata, key=shapes_key)
        geo_df = to_intrinsic(sdata, geo_df, points_key)

    gene_column = gene_column or get_feature_key(points, raise_error=True)

    log.info(f"Aggregating transcripts over {len(geo_df)} cells")
    return _count_transcripts_aligned(geo_df, points, gene_column)


def _count_transcripts_aligned(geo_df: gpd.GeoDataFrame, points: dd.DataFrame, value_key: str) -> AnnData:
    """Count transcripts per cell. The cells and points have to be aligned (i.e., in the same coordinate system)

    Args:
        geo_df: Cells geometries
        points: Transcripts dataframe
        value_key: Key of `points` containing the genes names

    Returns:
        An `AnnData` object of shape `(n_cells, n_genes)` with the counts per cell
    """
    points[value_key] = points[value_key].astype("category").cat.as_known()
    gene_names = points[value_key].cat.categories.astype(str)

    X = coo_matrix((len(geo_df), len(gene_names)), dtype=int)
    adata = AnnData(X=X, var=pd.DataFrame(index=gene_names))
    adata.obs_names = geo_df.index.astype(str)

    geo_df = geo_df.reset_index()

    X_partitions = []

    with ProgressBar():
        points.map_partitions(
            partial(_add_coo, X_partitions, geo_df, gene_column=value_key, gene_names=gene_names),
            meta=(),
        ).compute()

    for X_partition in X_partitions:
        adata.X += X_partition

    adata.X = adata.X.tocsr()

    if settings.gene_exclude_pattern is not None:
        adata = adata[:, ~adata.var_names.str.match(settings.gene_exclude_pattern, case=False, na=False)].copy()

    return adata


def _add_coo(
    X_partitions: list[coo_matrix],
    geo_df: gpd.GeoDataFrame,
    partition: pd.DataFrame,
    gene_column: str,
    gene_names: list[str],
) -> None:
    if settings.gene_exclude_pattern is not None:
        partition = partition[~partition[gene_column].str.match(settings.gene_exclude_pattern, case=False, na=False)]

    points_gdf = gpd.GeoDataFrame(partition, geometry=gpd.points_from_xy(partition["x"], partition["y"]))
    joined = geo_df.sjoin(points_gdf)
    cells_indices, column_indices = joined.index, joined[gene_column].cat.codes

    cells_indices = cells_indices[column_indices >= 0]
    column_indices = column_indices[column_indices >= 0]

    X_partition = coo_matrix(
        (np.full(len(cells_indices), 1), (cells_indices, column_indices)),
        shape=(len(geo_df), len(gene_names)),
    )

    X_partitions.append(X_partition)
