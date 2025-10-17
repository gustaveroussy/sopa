import logging

import geopandas as gpd
import numpy as np
from anndata import AnnData
from scipy.sparse import coo_matrix, csr_matrix
from spatialdata import SpatialData
from tqdm import tqdm

from ..shapes import expand_radius
from ..utils import to_intrinsic

log = logging.getLogger(__name__)


def aggregate_bins(
    sdata: SpatialData,
    shapes_key: str,
    bins_key: str,
    expand_radius_ratio: float = 0,
    no_overlap: bool = False,
) -> AnnData:
    """Aggregate bins (for instance, from Visium HD data) into cells.

    Args:
        sdata: The `SpatialData` object
        shapes_key: Key of the shapes containing the cell boundaries
        bins_key: Key of the table containing the bin-by-gene counts
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate bins from the cytoplasm.
        no_overlap: If `True`, bins belonging to multiple cells will be assigned to only one, based on transcript-profile proximity.

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

    indices_matrix = csr_matrix(
        (np.full(len(bin_within_cell), 1), (bin_within_cell["index_right"], bin_within_cell.index)),
        shape=(len(cells), len(bins)),
    )

    if no_overlap:
        log.warning("Unique bin assignments is currently experimental. Any feedback on GitHub is welcome.")
        indices_matrix = _get_unique_bins_assignments(indices_matrix, bins_table)

    adata = AnnData(indices_matrix @ bins_table.X, obs=cells[[]], var=bins_table.var)
    adata.obsm["spatial"] = np.stack([cells.centroid.x, cells.centroid.y], axis=1)
    adata.obsm["bins_assignments"] = indices_matrix
    return adata


def _get_unique_bins_assignments(
    indices_matrix: csr_matrix,
    bins_table: AnnData,
    n_components: int = 50,
) -> csr_matrix:
    shared_bins: np.ndarray = (indices_matrix.sum(0) >= 2).A1

    if not shared_bins.any():
        return indices_matrix

    indices_matrix: coo_matrix = indices_matrix.tocoo()
    no_overlap: coo_matrix = indices_matrix.copy()
    no_overlap.data[shared_bins[no_overlap.col]] = 0
    no_overlap.eliminate_zeros()
    no_overlap = no_overlap.tolil()

    adata_unique_map = AnnData(no_overlap @ bins_table.X)

    X_unique, X_shared = _pca_representation(n_components, adata_unique_map.X, bins_table.X[shared_bins])

    for bin_index_among_shared, bin_index in enumerate(tqdm(np.where(shared_bins)[0])):
        cell_indices = indices_matrix.row[indices_matrix.col == bin_index]

        distances: np.ndarray = ((X_unique[cell_indices] - X_shared[bin_index_among_shared]) ** 2).sum(1)
        cell_index = cell_indices[distances.argmin()]

        no_overlap[cell_index, bin_index] = 1

    return no_overlap.tocsr()


def _pca_representation(
    n_components: int, adata_unique_map: np.ndarray | csr_matrix, shared_bins: np.ndarray | csr_matrix
) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.decomposition import PCA

    n_components = min(n_components, shared_bins.shape[1] - 1, shared_bins.shape[0] - 1, adata_unique_map.shape[0] - 1)
    pca = PCA(n_components=n_components)

    X_unique = pca.fit_transform(_log1p_tpm(adata_unique_map))
    X_shared = pca.transform(_log1p_tpm(shared_bins))

    return X_unique, X_shared


def _log1p_tpm(x: np.ndarray | csr_matrix) -> np.ndarray:
    if isinstance(x, np.ndarray):
        x = x / (x.sum(axis=1, keepdims=True) + 1e-8) * 1e6
    else:
        x = x / (x.sum(axis=1) + 1e-8) * 1e6
    return np.log1p(x)
