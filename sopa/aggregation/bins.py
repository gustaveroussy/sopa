import logging

import geopandas as gpd
import numpy as np
from anndata import AnnData
from scipy.sparse import coo_matrix, lil_matrix
from spatialdata import SpatialData

from ..segmentation.shapes import expand_radius
from ..utils import to_intrinsic

log = logging.getLogger(__name__)


def aggregate_bins(
    sdata: SpatialData,
    shapes_key: str,
    bins_key: str,
    expand_radius_ratio: float = 0,
    unique_mapping: bool = False,
) -> AnnData:
    """Aggregate bins (for instance, from Visium HD data) into cells.

    Args:
        sdata: The `SpatialData` object
        shapes_key: Key of the shapes containing the cell boundaries
        bins_key: Key of the table containing the bin-by-gene counts
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate bins from the cytoplasm.
        unique_mapping: If `True`, bins belonging to multiples cells with be assigned to only one, based on transcript-profile proximity.

    Returns:
        An `AnnData` object of shape with the cell-by-gene count matrix
    """
    bins_table: AnnData = sdata.tables[bins_key]

    bins_shapes_key = sdata.get_annotated_regions(bins_table)
    bins_shapes_key = bins_shapes_key[0] if isinstance(bins_shapes_key, list) else bins_shapes_key
    bins = sdata.shapes[bins_shapes_key].loc[sdata.get_instance_key_column(bins_table).values]
    bins = gpd.GeoDataFrame(geometry=bins.centroid.values)  # bins as points

    cells = to_intrinsic(sdata, shapes_key, bins_shapes_key).reset_index(drop=True)
    cells = expand_radius(cells, expand_radius_ratio, no_overlap=False)

    bin_within_cell = gpd.sjoin(bins, cells)

    indices_matrix = coo_matrix(
        (np.full(len(bin_within_cell), 1), (bin_within_cell["index_right"], bin_within_cell.index)),
        shape=(len(cells), len(bins)),
    )

    if unique_mapping:
        indices_matrix = _get_unique_bins_mapping(indices_matrix, bins_table)

    adata = AnnData(indices_matrix @ bins_table.X, obs=cells[[]], var=bins_table.var)
    adata.obsm["spatial"] = np.stack([cells.centroid.x, cells.centroid.y], axis=1)
    return adata


def _get_unique_bins_mapping(
    indices_matrix: coo_matrix,
    bins_table: AnnData,
    n_components: int = 50,
) -> coo_matrix | lil_matrix:
    shared_bins: np.ndarray = (indices_matrix.sum(0) >= 2).A1

    if not shared_bins.any():
        return indices_matrix

    unique_mapping: coo_matrix = indices_matrix.copy()
    unique_mapping.data[shared_bins[unique_mapping.col]] = 0
    unique_mapping.eliminate_zeros()
    unique_mapping = unique_mapping.tolil()

    adata_unique_mapping = AnnData(unique_mapping @ bins_table.X)

    X_unique, X_shared = _pca_representation(n_components, adata_unique_mapping.X, bins_table.X[shared_bins])

    for bin_index_among_shared, bin_index in enumerate(np.where(shared_bins)[0]):
        cell_indices = indices_matrix.row[indices_matrix.col == bin_index]

        distances: np.ndarray = ((X_unique[cell_indices] - X_shared[bin_index_among_shared]) ** 2).sum(1)
        cell_index = cell_indices[distances.argmin()]

        unique_mapping[cell_index, bin_index] = 1

    return unique_mapping


def _pca_representation(
    n_components: int, table_unique_mapping: np.ndarray, shared_table: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.decomposition import PCA

    n_components = min(n_components, shared_table.shape[1] - 1)
    pca = PCA(n_components=n_components)

    X_unique = pca.fit_transform(_log1p_tpm(table_unique_mapping))
    X_shared = pca.transform(_log1p_tpm(shared_table))

    return X_unique, X_shared


def _log1p_tpm(x: np.ndarray) -> np.ndarray:
    x = x / (x.sum(axis=1, keepdims=True) * 1e6 + 1e-6)
    return np.log1p(x)
