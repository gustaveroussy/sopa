"""
Copied from squidpy (to not have squidpy as a dependency)
Functions for building graphs from spatial coordinates.
"""

from __future__ import annotations

import logging
import warnings
from functools import partial
from itertools import chain
from typing import Iterable

import numpy as np
from anndata import AnnData
from anndata.utils import make_index_unique
from scipy.sparse import SparseEfficiencyWarning, block_diag, csr_matrix, spmatrix
from scipy.spatial import Delaunay
from sklearn.metrics.pairwise import euclidean_distances
from spatialdata import SpatialData

from .._constants import SopaKeys

log = logging.getLogger(__name__)
__all__ = ["spatial_neighbors"]


def spatial_neighbors(
    adata: AnnData | SpatialData,
    radius: tuple[float, float] | None,
    library_key: str | None = None,
    percentile: float | None = None,
    set_diag: bool = False,
):
    """Create a Delaunay graph from spatial coordinates. This function comes from [squidpy](https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.spatial_neighbors.html#squidpy.gr.spatial_neighbors).

    Args:
        adata: AnnData object
        radius: tuple that prunes the final graph to only contain edges in interval `[min(radius), max(radius)]`. If `None`, all edges are kept.
        library_key: Optional batch key in adata.obs
        percentile: Percentile of the distances to use as threshold.
        set_diag: Whether to set the diagonal of the spatial connectivities to `1.0`.
    """
    if isinstance(adata, SpatialData):
        adata = adata.tables[SopaKeys.TABLE]

    assert (
        radius is None or len(radius) == 2
    ), "Radius is expected to be a tuple (min_radius, max_radius)"

    log.info("Computing delaunay graph")

    if library_key is not None:
        assert adata.obs[library_key].dtype == "category"
        libs = adata.obs[library_key].cat.categories
        make_index_unique(adata.obs_names)
    else:
        libs = [None]

    _build_fun = partial(
        _spatial_neighbor,
        set_diag=set_diag,
        radius=radius,
        percentile=percentile,
    )

    if library_key is not None:
        mats: list[tuple[spmatrix, spmatrix]] = []
        ixs = []  # type: ignore[var-annotated]
        for lib in libs:
            ixs.extend(np.where(adata.obs[library_key] == lib)[0])
            mats.append(_build_fun(adata[adata.obs[library_key] == lib]))
        ixs = np.argsort(ixs)  # type: ignore[assignment] # invert
        Adj = block_diag([m[0] for m in mats], format="csr")[ixs, :][:, ixs]
        Dst = block_diag([m[1] for m in mats], format="csr")[ixs, :][:, ixs]
    else:
        Adj, Dst = _build_fun(adata)

    adata.obsp["spatial_connectivities"] = Adj
    adata.obsp["spatial_distances"] = Dst


def _spatial_neighbor(
    adata: AnnData,
    radius: float | tuple[float, float] | None = None,
    set_diag: bool = False,
    percentile: float | None = None,
) -> tuple[csr_matrix, csr_matrix]:
    coords = adata.obsm["spatial"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SparseEfficiencyWarning)
        Adj, Dst = _build_connectivity(
            coords,
            set_diag=set_diag,
        )

    if isinstance(radius, Iterable):
        minn, maxx = sorted(radius)[:2]  # type: ignore[var-annotated]
        mask = (Dst.data < minn) | (Dst.data > maxx)
        a_diag = Adj.diagonal()

        Dst.data[mask] = 0.0
        Adj.data[mask] = 0.0
        Adj.setdiag(a_diag)

    if percentile is not None:
        threshold = np.percentile(Dst.data, percentile)
        Adj[Dst > threshold] = 0.0
        Dst[Dst > threshold] = 0.0

    Adj.eliminate_zeros()
    Dst.eliminate_zeros()

    return Adj, Dst


def _build_connectivity(
    coords: np.ndarray,
    set_diag: bool = False,
) -> csr_matrix | tuple[csr_matrix, csr_matrix]:
    N = coords.shape[0]

    tri = Delaunay(coords)
    indptr, indices = tri.vertex_neighbor_vertices
    Adj = csr_matrix((np.ones_like(indices, dtype=np.float64), indices, indptr), shape=(N, N))

    # fmt: off
    dists = np.array(list(chain(*(
        euclidean_distances(coords[indices[indptr[i] : indptr[i + 1]], :], coords[np.newaxis, i, :])
        for i in range(N)
        if len(indices[indptr[i] : indptr[i + 1]])
    )))).squeeze()
    Dst = csr_matrix((dists, indices, indptr), shape=(N, N))
    # fmt: on

    # radius-based filtering needs same indices/indptr: do not remove 0s
    Adj.setdiag(1.0 if set_diag else Adj.diagonal())
    Dst.setdiag(0.0)

    return Adj, Dst


def _check_has_delaunay(adata: AnnData):
    message = " key not in adata.obsp, consider running the delaunay graph (e.g., `from sopa.spatial import spatial_neighbors; spatial_neighbors(adata, [0, 40])`)"
    assert "spatial_connectivities" in adata.obsp, "spatial_connectivities" + message
    assert "spatial_distances" in adata.obsp, "spatial_distances" + message
