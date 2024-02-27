from __future__ import annotations


class FileNames:
    IMAGE = "morphology.ome.tif"
    POINTS = "transcripts.zarr.zip"
    SHAPES = "cells.zarr.zip"
    H5AD = "adata.h5ad"
    TABLE = "cell_feature_matrix.zarr.zip"
    CELL_CATEGORIES = "analysis.zarr.zip"
    METADATA = "experiment.xenium"


class ExplorerConstants:
    GRID_SIZE = 250
    TILE_SIZE = 1024
    QUALITY_SCORE = 40
    MICRONS_TO_PIXELS = 4.705882
    PIXELS_TO_MICRONS = 0.2125

    COLORS = ["white", 400, 500, 600, 700]
    NUCLEUS_COLOR = "white"
    KNOWN_CHANNELS = {"DAPI": "white", "DNA1": "white", "DNA2": "white", "DAPI (000)": "white"}


class Versions:
    EXPERIMENT = [2, 0]
    GROUPS = [5, 0]
    CELL_CATEGORIES = [1, 0]


def cell_categories_attrs() -> dict:
    return {
        "major_version": Versions.CELL_CATEGORIES[0],
        "minor_version": Versions.CELL_CATEGORIES[1],
        "number_groupings": 0,
        "grouping_names": [],
        "group_names": [],
    }


def cell_summary_attrs() -> dict:
    return {
        "column_descriptions": [
            "Cell centroid in X",
            "Cell centroid in Y",
            "Cell area",
            "Nucleus centroid in X",
            "Nucleus centroid in Y",
            "Nucleus area",
            "z_level",
        ],
        "column_names": [
            "cell_centroid_x",
            "cell_centroid_y",
            "cell_area",
            "nucleus_centroid_x",
            "nucleus_centroid_y",
            "nucleus_area",
            "z_level",
        ],
    }


def group_attrs() -> dict:
    return {
        "major_version": Versions.GROUPS[0],
        "minor_version": Versions.GROUPS[1],
        "name": "CellSegmentationDataset",
        "polygon_set_descriptions": [
            "NA",
            "NA",
        ],
        "polygon_set_display_names": ["Nucleus boundaries", "Cell boundaries"],
        "polygon_set_names": ["nucleus", "cell"],
        "spatial_units": "microns",
    }


def experiment_dict(run_name: str, region_name: str, num_cells: int, pixel_size: float) -> dict:
    return {
        "major_version": Versions.EXPERIMENT[0],
        "minor_version": Versions.EXPERIMENT[1],
        "run_name": run_name,
        "region_name": region_name,
        "experiment_uuid": "N/A",
        "panel_tissue_type": "N/A",
        "run_start_time": "N/A",
        "preservation_method": "N/A",
        "num_cells": num_cells,
        "transcripts_per_cell": 0,
        "transcripts_per_100um": 0,
        "cassette_name": "N/A",
        "slide_id": "N/A",
        "panel_design_id": "N/A",
        "panel_name": "N/A",
        "panel_organism": "Human",
        "panel_num_targets_predesigned": 0,
        "panel_num_targets_custom": 0,
        "pixel_size": pixel_size,
        "instrument_sn": "N/A",
        "instrument_sw_version": "N/A",
        "analysis_sw_version": "xenium-1.3.0.5",
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


def image_metadata(channel_names: list[str], pixel_size: float) -> dict:
    return {
        "SignificantBits": 8,
        "PhysicalSizeX": pixel_size,
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": pixel_size,
        "PhysicalSizeYUnit": "µm",
        "Channel": {"Name": channel_names},
    }
