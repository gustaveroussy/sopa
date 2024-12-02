import json
import logging
import os
import sys
import warnings
from functools import partial
from pathlib import Path

import geopandas as gpd
import numpy as np
from spatialdata import SpatialData

from ... import settings
from ..._constants import SopaAttrs, SopaKeys
from ...utils import (
    delete_transcripts_patches_dirs,
    get_transcripts_patches_dirs,
    to_intrinsic,
)
from .._transcripts import _check_transcript_patches, resolve

log = logging.getLogger(__name__)


def comseg(
    sdata: SpatialData,
    config: dict | str | None = None,
    min_area: float = 0,
    delete_cache: bool = True,
    recover: bool = False,
    key_added: str = SopaKeys.COMSEG_BOUNDARIES,
    patch_index: int | None = None,
):
    """Run [ComSeg](https://comseg.readthedocs.io/en/latest/) segmentation on a SpatialData object, and add a GeoDataFrame containing the cell boundaries.

    !!! warning "ComSeg installation"
        Make sure to install ComSeg (`pip install comseg`) for this method to work.

    !!! info "Transcript patches"
        To use ComSeg, make sure to run [sopa.make_transcript_patches][] with a `prior_shapes_key` and `write_cells_centroids=True`.

    Args:
        sdata: A `SpatialData` object.
        config: Optional configuration dictionary or path to a JSON file containing a valid ComSeg config. By default, a configuration is inferred based on the cell area of the prior segmentation.
        min_area: Minimal area (in microns^2) of a cell to be considered.
        delete_cache: Whether to delete the cache after segmentation.
        recover: If `True`, recover the cache from a failed segmentation, and continue.
        key_added: Name of the shapes element to be added to `sdata`.
        patch_index: Index of the patch to segment (we do not recommend to set this argument). By default, segment all patches.
    """
    _check_transcript_patches(sdata, with_prior=True)

    if config is None or not len(config):
        config = _get_default_config(sdata, sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES])
    elif isinstance(config, str):
        with open(config, "r") as f:
            config = json.load(f)

    assert "gene_column" in config, "'gene_column' not found in config"

    config["prior_name"] = sdata[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.PRIOR_SHAPES_KEY].iloc[0]

    if patch_index is not None:
        patch_dir = Path(sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES].loc[patch_index, SopaKeys.CACHE_PATH_KEY])
        comseg_patch(patch_dir, config, recover)
        return

    patches_dirs = get_transcripts_patches_dirs(sdata)

    _functions = [partial(comseg_patch, patch_dir, config, recover) for patch_dir in patches_dirs]
    settings._run_with_backend(_functions)

    resolve(sdata, patches_dirs, config["gene_column"], min_area=min_area, key_added=key_added)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added

    if delete_cache:
        delete_transcripts_patches_dirs(sdata)


def comseg_patch(patch_dir: Path, config: dict, recover: bool = False):
    import json

    try:
        import comseg
        from comseg import dataset as ds
        from comseg import dictionary
    except ModuleNotFoundError:
        raise ModuleNotFoundError("Install comseg (`pip install comseg`) for this method to work")

    assert comseg.__version__ >= "1.8.2", "comseg version should be >= 1.8.2"

    if (
        recover
        and (patch_dir / "segmentation_counts.h5ad").exists()
        and (patch_dir / "segmentation_polygons.json").exists()
    ):
        return

    if "disable_tqdm" not in config:
        config["disable_tqdm"] = True

    warnings.filterwarnings("ignore", message="param_sctransform is none")
    warnings.filterwarnings("ignore", message="Series.__getitem__")

    with HiddenPrints():
        dataset = ds.ComSegDataset(
            path_dataset_folder=patch_dir,
            dict_scale=config["dict_scale"],
            mean_cell_diameter=config["mean_cell_diameter"],
            gene_column=config["gene_column"],
            image_csv_files=["transcripts.csv"],
            centroid_csv_files=["centroids.csv"],
            path_cell_centroid=patch_dir,
            min_nb_rna_patch=config.get("min_nb_rna_patch", 0),
            prior_name=config["prior_name"],
            disable_tqdm=config["disable_tqdm"],
        )

        dataset.compute_edge_weight(config=config)

        Comsegdict = dictionary.ComSegDict(
            dataset=dataset,
            mean_cell_diameter=config["mean_cell_diameter"],
            disable_tqdm=config["disable_tqdm"],
        )

        Comsegdict.run_all(config=config)

        assert config.get("return_polygon", True), "Only return_polygon=True is supported in sopa"

        anndata_comseg, json_dict = Comsegdict.anndata_from_comseg_result(config=config)

        anndata_comseg.write_h5ad(patch_dir / "segmentation_counts.h5ad")
        with open(patch_dir / "segmentation_polygons.json", "w") as f:
            json.dump(json_dict["transcripts"], f)


def _get_default_config(sdata: SpatialData, patches: gpd.GeoDataFrame) -> dict:
    prior_shapes_key = patches[SopaKeys.PRIOR_SHAPES_KEY].iloc[0]
    points_key = patches[SopaKeys.POINTS_KEY].iloc[0]

    cells_area = to_intrinsic(sdata, prior_shapes_key, points_key).area
    mean_cell_diameter = 2 * np.sqrt(cells_area / np.pi).mean()

    config = {
        "dict_scale": {"x": 1, "y": 1, "z": 1},  # spot coordinates already in µm
        "mean_cell_diameter": mean_cell_diameter,
        "max_cell_radius": mean_cell_diameter * 1.75,
        "norm_vector": False,
        "alpha": 0.5,  # alpha value to compute the polygon https://pypi.org/project/alphashape/
        "allow_disconnected_polygon": False,
        "min_rna_per_cell": 20,  # minimal number of RNAs for a cell to be taken into account
        "gene_column": "genes",
    }

    log.info(f"The Comseg config was not provided, using the following by default:\n{config}")

    return config


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
