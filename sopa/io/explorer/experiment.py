import json

EXPERIMENT = {
    "major_version": 2,
    "minor_version": 0,
    "run_start_time": "N/A",
    "preservation_method": "N/A",
    "num_cells": 1,
    "transcripts_per_cell": 0,
    "transcripts_per_100um": 0,
    "cassette_name": "N/A",
    "slide_id": "N/A",
    "panel_design_id": "N/A",
    "panel_name": "R&D Panel",
    "panel_organism": "Human",
    "panel_num_targets_predesigned": 254,
    "panel_num_targets_custom": 65,
    "pixel_size": 0.2125,
    "instrument_sn": "R&D",
    "instrument_sw_version": "R&D",
    "analysis_sw_version": "xenium-1.3.0.5",
    "experiment_uuid": "",
    "cassette_uuid": "",
    "roi_uuid": "",
    "z_step_size": 3.0,
    "well_uuid": "",
    "calibration_uuid": "N/A",
    "images": {
        "morphology_filepath": "morphology.ome.tif",
        "morphology_mip_filepath": "morphology_mip.ome.tif",
        "morphology_focus_filepath": "morphology_focus.ome.tif",
    },
    "xenium_explorer_files": {
        "transcripts_zarr_filepath": "transcripts.zarr.zip",
        "cells_zarr_filepath": "cells.zarr.zip",
        "cell_features_zarr_filepath": "cell_feature_matrix.zarr.zip",
        "analysis_zarr_filepath": "analysis.zarr.zip",
        "analysis_summary_filepath": "analysis_summary.html",
    },
}


def write_experiment(path, run_name, tissue, region_name, uuid):
    EXPERIMENT["run_name"] = run_name
    EXPERIMENT["panel_tissue_type"] = tissue
    EXPERIMENT["region_name"] = region_name
    EXPERIMENT["experiment_uuid"] = uuid

    with open(path, "w") as f:
        json.dump(EXPERIMENT, f, indent=4)
