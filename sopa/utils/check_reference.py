import logging

import anndata
import click
from scipy.sparse import issparse

MIN_GENES = 100
MIN_CELLS = 5_000
CT_KEY = "ct_level"
log = logging.getLogger(__name__)

# TODO: improve and remove click


@click.command()
@click.option("-p", "--path", type=str, help="Path to the h5ad reference")
def main(path: str):
    adata = anndata.read_h5ad(path)

    assert (
        adata.n_vars >= MIN_GENES
    ), f"The reference must have at least {MIN_GENES} genes. Found {adata.n_vars}."

    assert (
        adata.n_obs >= MIN_CELLS
    ), f"The reference must have at least {MIN_CELLS} cells. Found {adata.n_obs}."

    assert (
        "counts" in adata.layers
    ), "No 'counts' found in adata.layers. This must be a sparse matrix."

    assert issparse(adata.layers["counts"]), "The 'counts' layer must be sparse."

    assert (
        f"{CT_KEY}0" in adata.obs
    ), f"Provide at least one level of cell-type annotation, inside adata.obs['{CT_KEY}0']"

    max_level = 1
    while f"{CT_KEY}{max_level}" in adata.obs:
        ct_parent, ct_children = f"{CT_KEY}{max_level - 1}", f"{CT_KEY}{max_level}"
        counts = adata.obs.groupby(ct_children, observed=True)[ct_parent].value_counts().unstack()
        n_parents = (counts > 0).sum(1)
        assert (
            n_parents == 1
        ).all(), f"All populations on {ct_children} must have one and only one parent in {ct_parent}. The number of parents is the following:\n{n_parents}"
        max_level += 1

    log.info(f"{max_level} level(s) of annotations were found")

    max_x = adata.X.max()
    assert (
        max_x != int(max_x) and max_x < 30
    ), "Please preprocess the expressions with sc.pp.normalize_total and sc.pp.log1p"

    log.info("The reference passed the sanity check!")


if __name__ == "__main__":
    main()
