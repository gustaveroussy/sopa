import typer

app_check = typer.Typer()
option = typer.Option()


def _open_config(path: str) -> dict:
    import logging

    import yaml

    log = logging.getLogger(__name__)

    try:
        with open(path, "r") as f:
            return yaml.safe_load(f)
    except:
        log.warn(
            f"Config file '{path}' could't be read. Make sure that the file exist and that it is a YAML file"
        )
        return


@app_check.command()
def reference(reference_path: str, cell_type_key: str = "cell_type"):
    """Perform sanity checks on a tangram scRNAseq reference

    Args:\n
        reference_path: Path to the scRNAseq reference (usually, a h5ad file)\n
        cell_type_key: Key of adata.obs containing the cell types
    """
    import logging

    import anndata

    log = logging.getLogger(__name__)

    adata = anndata.read(reference_path)

    MIN_GENES = 100
    MIN_CELLS = 1_000

    if adata.n_vars >= MIN_GENES:
        log.warn(f"The reference must have at least {MIN_GENES} genes. Found {adata.n_vars}.")

    if adata.n_obs >= MIN_CELLS:
        log.warn(f"The reference must have at least {MIN_CELLS} cells. Found {adata.n_obs}.")

    if not (cell_type_key in adata.obs):
        log.warn(
            f"Column adata.obs['{cell_type_key}'] not found. Update your anndata object, or provide another --cell-type-key argument"
        )

    level = 1
    current_key = cell_type_key
    next_key = f"{cell_type_key}_level{level}"
    while next_key in adata.obs:
        counts = adata.obs.groupby(next_key, observed=True)[current_key].value_counts().unstack()
        n_parents = (counts > 0).sum(1)
        if (n_parents == 1).all():
            log.warn(
                f"All populations on {next_key} must have one and only one parent in {current_key}. The number of parents is the following:\n{n_parents}"
            )

        level += 1
        current_key, next_key = next_key, f"{cell_type_key}_level{level}"

    log.info(
        f"{level} level(s) of annotations were found. To provide more annotation levels, add adata.obs['{cell_type_key}_level{level + 1}']"
    )

    log.info("Reference sanity check completed")


@app_check.command()
def config(path: str):
    """Perform sanity checks on a sopa yaml config

    Args:\n
        path: Path to the YAML config\n
    """
    import logging

    log = logging.getLogger(__name__)

    config = _open_config(path)

    REQUIRED_KEYS = ["read", "patchify", "segmentation", "aggregate"]

    for key in REQUIRED_KEYS:
        if not key in config:
            log.warn(f"Key '{key}' is missing from the config")
        elif not isinstance(config[key], dict):
            log.warn(
                f"config['{key}'] must be a dictionary (i.e. YAML indentation). But found {type(config[key])}"
            )

    # TODO: more complete check

    log.info("Config sanity check completed")
