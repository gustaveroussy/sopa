import anndata


def to_adata(path_loom: str):
    adata = anndata.read_loom(path_loom, obs_names="Name", var_names="Name")
