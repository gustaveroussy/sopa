from __future__ import annotations

from typing import Callable

import anndata
import geopandas as gpd
import numpy as np
import scanpy as sc
from spatial_image import SpatialImage
from spatialdata import SpatialData

from sopa._constants import SopaKeys


def leiden_clustering(X: np.ndarray, **kwargs):
    adata = anndata.AnnData(X=X)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, **kwargs)
    return adata.obs["leiden"].values


METHODS_DICT = {
    "leiden": leiden_clustering,
}


def cluster_embeddings(
    sdata: SpatialData,
    element: SpatialImage | str,
    method: Callable | str = "leiden",
    key_added: str = "cluster",
    **method_kwargs: str,
) -> gpd.GeoDataFrame:
    """Cluster the patches embeddings using a clustering method

    Args:
        sdata: A `SpatialData` object
        element: The `SpatialImage` containing the embeddings, or the name of the element
        method: Callable that takes as an input an array of size `(n_patches x embedding_size)` and returns an array of clusters of size `n_patches`, or an available method name (`leiden`)
        key_added: The key containing the clusters to be added to the patches `GeoDataFrame`
        method_kwargs: kwargs provided to the method callable

    Returns:
        The patches `GeoDataFrame` with a new column `key_added` containing the patches clusters
    """
    if isinstance(element, str):
        element = sdata.images[element]

    if isinstance(method, str):
        assert (
            method in METHODS_DICT
        ), f"Method {method} is not available. Use one of: {', '.join(METHODS_DICT.keys())}"
        method = METHODS_DICT[method]

    gdf_patches = sdata[SopaKeys.EMBEDDINGS_PATCHES_KEY]

    ilocs = np.array(list(gdf_patches.ilocs))
    embeddings = element.compute().data[:, ilocs[:, 1], ilocs[:, 0]].T

    gdf_patches[key_added] = method(embeddings, **method_kwargs)

    return gdf_patches
