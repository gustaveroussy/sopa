from pathlib import Path

import loompy
from anndata import AnnData


def read_loom(path: Path) -> AnnData:
    ds = loompy.connect(path)
    adata: AnnData = AnnData(ds[:, :], dtype="int64").T
    adata.var_names = ds.ra.Name
    adata.obs_names = ds.ca.Name
    return adata
