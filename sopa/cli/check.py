from __future__ import annotations

import typer

app_check = typer.Typer()


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
def reference(
    reference_path: str = typer.Argument(help="Path to the scRNAseq reference as a `.h5ad` file"),
    cell_type_key: str = typer.Option(help="Key of adata.obs containing the cell types"),
):
    """Perform sanity checks on a tangram scRNAseq reference"""
    import logging

    import anndata

    log = logging.getLogger(__name__)

    try:
        adata = anndata.read_h5ad(reference_path)
    except:
        log.warn(
            f"scRNAseq reference at '{reference_path}' could't be read. Make sure that the file exist and that it can be open by `anndata.read_h5ad`"
        )
        return

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


def _get(config, *args):
    for arg in args:
        if arg not in config:
            return False
        config = config[arg]
    return config


def _check_dict(config: dict, d: dict, log, prefix: str = "config"):
    for key, values in d.items():
        if key not in config:
            log.warn(f"Required config key {prefix}['{key}'] not found")
        elif isinstance(values, dict):
            _check_dict(config[key], values, log, f"{prefix}['{key}']")
        elif isinstance(values, list):
            for element in values:
                if isinstance(element, str) and element in config[key]:
                    break
                if all(x in config[key] for x in element):
                    break
            else:
                display = "\n  - ".join(
                    element if isinstance(element, str) else " AND ".join(element)
                    for element in values
                )
                log.warn(f"One of these element must be in {prefix}['{key}']:\n  - {display}")


CONFIG_REQUIREMENTS = {
    "read": ["technology"],
    "patchify": [
        ["patch_width_pixel", "patch_overlap_pixel"],
        ["patch_width_microns", "patch_overlap_microns"],
    ],
    "segmentation": ["cellpose", "baysor"],
}


@app_check.command()
def config(path: str = typer.Argument(help="Path to the YAML config")):
    """Perform sanity checks on a sopa yaml config"""
    import logging

    log = logging.getLogger(__name__)

    config = _open_config(path)

    _check_dict(config, CONFIG_REQUIREMENTS, log)

    if _get(config, "annotation", "method") == "tangram":
        sc_reference_path = _get(config, "annotation", "args", "sc_reference_path")
        if not sc_reference_path:
            log.warn("Tangram used but no config['annotation']['args']['sc_reference_path'] found")
        else:
            cell_type_key = _get(config, "annotation", "args", "cell_type_key") or "cell_type"
            reference(sc_reference_path, cell_type_key=cell_type_key)

    log.info("Minimal config sanity check completed")
