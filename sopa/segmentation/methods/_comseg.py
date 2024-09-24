from __future__ import annotations

from pathlib import Path

from ..._constants import SopaKeys


def comseg_patch(temp_dir: str, patch_index: int, config: dict):
    import json

    try:
        import comseg
        from comseg import dataset as ds
        from comseg import dictionary
    except ModuleNotFoundError:
        raise ModuleNotFoundError("Install comseg (`pip install comseg`) for this method to work")

    assert comseg.__version__ >= "1.3", "comseg version should be >= 1.3"

    path_dataset_folder = Path(temp_dir) / str(patch_index)

    dataset = ds.ComSegDataset(
        path_dataset_folder=path_dataset_folder,
        dict_scale=config["dict_scale"],
        mean_cell_diameter=config["mean_cell_diameter"],
        gene_column=config["gene_column"],
        image_csv_files=["transcripts.csv"],
        centroid_csv_files=["centroids.csv"],
        path_cell_centroid=path_dataset_folder,
        min_nb_rna_patch=config.get("min_nb_rna_patch", 0),
        prior_name=config.get("prior_name", SopaKeys.DEFAULT_CELL_KEY),
    )

    dataset.compute_edge_weight(config=config)

    Comsegdict = dictionary.ComSegDict(
        dataset=dataset,
        mean_cell_diameter=config["mean_cell_diameter"],
    )

    Comsegdict.run_all(config=config)

    if "return_polygon" in config:
        assert config["return_polygon"] is True, "Only return_polygon=True is supported in sopa"
    anndata_comseg, json_dict = Comsegdict.anndata_from_comseg_result(config=config)
    anndata_comseg.write_h5ad(path_dataset_folder / "segmentation_counts.h5ad")
    with open(path_dataset_folder / "segmentation_polygons.json", "w") as f:
        json.dump(json_dict["transcripts"], f)
