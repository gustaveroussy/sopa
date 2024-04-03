from __future__ import annotations

import ast

import typer

from .utils import SDATA_HELPER

app_annotate = typer.Typer()


@app_annotate.command()
def fluorescence(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    marker_cell_dict: str = typer.Option(callback=ast.literal_eval),
    cell_type_key: str = typer.Option(
        "cell_type", help="Key added in `adata.obs` corresponding to the cell type"
    ),
):
    """Simple annotation based on fluorescence, where each provided channel corresponds to one cell type.

    For each cell, one z-score statistic is computed and the population with the highest z-score is attributed.
    """
    from sopa._constants import SopaKeys
    from sopa._sdata import save_table
    from sopa.annotation.fluorescence import higher_z_score
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)

    assert SopaKeys.TABLE in sdata.tables, f"No '{SopaKeys.TABLE}' found in sdata.tables"

    higher_z_score(sdata.tables[SopaKeys.TABLE], marker_cell_dict, cell_type_key)
    save_table(sdata, SopaKeys.TABLE)


@app_annotate.command()
def tangram(
    sdata_path: str = typer.Argument(help=SDATA_HELPER),
    sc_reference_path: str = typer.Option(
        help="Path to the scRNAseq annotated reference, as a `.h5ad` file"
    ),
    cell_type_key: str = typer.Option(help="Key of `adata_ref.obs` containing the cell-types"),
    reference_preprocessing: str = typer.Option(
        None,
        help="Preprocessing method applied to the reference. Either None (raw counts), or `normalized` (sc.pp.normalize_total) or `log1p` (sc.pp.normalize_total and sc.pp.log1p)",
    ),
    bag_size: int = typer.Option(
        10_000,
        help="Number of cells in each bag of the spatial table. Low values will decrease the memory usage",
    ),
    max_obs_reference: int = typer.Option(
        10_000,
        help="Maximum samples to be considered in the reference for tangram. Low values will decrease the memory usage",
    ),
):
    """Tangram segmentation (i.e., uses an annotated scRNAseq reference to transfer cell-types)"""
    import anndata

    from sopa._constants import SopaKeys
    from sopa._sdata import save_table
    from sopa.annotation.tangram.run import tangram_annotate
    from sopa.io.standardize import read_zarr_standardized

    sdata = read_zarr_standardized(sdata_path)
    adata_sc = anndata.read_h5ad(sc_reference_path)

    tangram_annotate(
        sdata,
        adata_sc,
        cell_type_key,
        reference_preprocessing=reference_preprocessing,
        bag_size=bag_size,
        max_obs_reference=max_obs_reference,
    )
    save_table(sdata, SopaKeys.TABLE)
