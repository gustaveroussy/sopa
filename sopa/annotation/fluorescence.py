import logging

import numpy as np
import pandas as pd
from anndata import AnnData

from .._constants import SopaKeys

log = logging.getLogger(__name__)


def preprocess_fluo(adata: AnnData) -> pd.DataFrame:
    if SopaKeys.INTENSITIES_OBSM in adata.obsm:
        df = adata.obsm[SopaKeys.INTENSITIES_OBSM]
    else:
        df = adata.to_df()

    divider = 5 * np.quantile(adata.X, 0.2, axis=0)
    divider[divider == 0] = adata.X.max(axis=0)[divider == 0]

    scaled = np.arcsinh(df / (5 * np.quantile(df, 0.2, axis=0)))
    return (scaled - scaled.mean(0)) / scaled.std(0)


def higher_z_score(adata: AnnData, marker_cell_dict: dict, cell_type_key: str = "cell_type"):
    adata.obsm[SopaKeys.Z_SCORES] = preprocess_fluo(adata)

    markers, cell_types = list(marker_cell_dict.keys()), np.array(list(marker_cell_dict.values()))
    ct_indices = adata.obsm[SopaKeys.Z_SCORES][markers].values.argmax(1)

    adata.obs[cell_type_key] = cell_types[ct_indices]
    adata.uns[SopaKeys.UNS_KEY][SopaKeys.UNS_CELL_TYPES] = [cell_type_key]

    log.info(f"Annotation counts: {adata.obs[cell_type_key].value_counts()}")
