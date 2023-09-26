import logging

import numpy as np
from anndata import AnnData

log = logging.getLogger(__name__)


def higher_z_score(adata: AnnData, marker_cell_dict: dict, key: str = "cell_type"):
    scaled = np.arcsinh(adata.X / (5 * np.quantile(adata.X, 0.2, axis=0)))

    adata.layers["z_scores"] = (scaled - scaled.mean(0)) / scaled.std(0)

    markers, cell_types = list(marker_cell_dict.keys()), np.array(list(marker_cell_dict.values()))
    ct_indices = adata[:, markers].layers["z_scores"].argmax(1)

    adata.obs[key] = cell_types[ct_indices]

    log.info(f"Annotation counts: {adata.obs[key].value_counts()}")
