from __future__ import annotations

from pathlib import Path

from .constants import STAINING_BASED_METHODS, TOY_READERS, TRANSCRIPT_BASED_METHODS


def validate_config(config: dict):
    """Basic config sanity check before running Snakemake"""
    assert "read" in config and "technology" in config["read"], (
        "Invalid config. Provide a 'read' section in the config file, and specify the 'technology' key"
    )

    check_segmentation_methods(config)

    check_data_paths(config)

    check_prior_shapes_key(config)

    if "baysor" in config["segmentation"]:
        check_baysor_executable_path(config)

    if "proseg" in config["segmentation"]:
        assert "prior_shapes_key" in config["segmentation"]["proseg"], (
            "Invalid config. Provide a 'prior_shapes_key' key in the 'proseg' section of the config file"
        )

        if "patch_width_microns" in config["patchify"]:
            assert config["patchify"]["patch_width_microns"] == -1, (
                "Invalid config. 'patch_width_microns' must be -1 for 'proseg' segmentation method"
            )
        config["patchify"]["patch_width_microns"] = -1

    return config


def check_segmentation_methods(config: dict):
    assert "segmentation" in config, "Invalid config. Provide a 'segmentation' section in the config file"

    assert len(config["segmentation"]) > 0, (
        "Invalid config. Provide at least one segmentation method in the 'segmentation' section"
    )

    assert sum(method in config["segmentation"] for method in TRANSCRIPT_BASED_METHODS) <= 1, (
        f"Only one of the following methods can be used: {TRANSCRIPT_BASED_METHODS}"
    )

    assert sum(method in config["segmentation"] for method in STAINING_BASED_METHODS) <= 1, (
        f"Only one of the following methods can be used: {STAINING_BASED_METHODS}"
    )

    if "stardist" in config["segmentation"]:
        assert not any(method in config["segmentation"] for method in TRANSCRIPT_BASED_METHODS), (
            "Invalid config. 'stardist' cannot be combined with transcript-based methods"
        )


def check_prior_shapes_key(config: dict):
    for method in TRANSCRIPT_BASED_METHODS:
        if method in config["segmentation"]:
            if "cellpose" in config["segmentation"]:
                config["segmentation"][method]["prior_shapes_key"] = "cellpose_boundaries"
            elif "cell_key" in config["segmentation"][method]:  # backward compatibility
                print("Snakemake argument 'cell_key' is deprecated. Use 'prior_shapes_key' instead.")
                config["segmentation"][method]["prior_shapes_key"] = config["segmentation"][method]["cell_key"]


def check_data_paths(config: dict):
    config["is_toy_reader"] = config["read"]["technology"] in TOY_READERS

    if config["is_toy_reader"]:
        config["data_path"] = "."

    assert "data_path" in config or "sdata_path" in config, (
        "Invalid config. Provide '--config data_path=...' when running the pipeline"
    )

    if "data_path" in config and "sdata_path" not in config:
        config["sdata_path"] = Path(config["data_path"]).with_suffix(".zarr")
        print(
            f"SpatialData object path set to default: {config['sdata_path']}\nTo change this behavior, provide `--config sdata_path=...` when running the snakemake pipeline"
        )

    if "data_path" not in config:
        assert Path(config["sdata_path"]).exists(), (
            f"When `data_path` is not provided, the spatial data object must exist, but the directory doesn't exists: {config['sdata_path']}"
        )
        config["data_path"] = []


def check_baysor_executable_path(config: dict) -> str:
    import shutil

    if shutil.which("baysor") is not None:
        return

    default_path = Path.home() / ".julia" / "bin" / "baysor"
    if default_path.exists():
        return

    if "executables" in config and "baysor" in config["executables"]:
        raise ValueError(
            "The config['executables']['baysor'] argument is deprecated. Please set a 'baysor' alias instead."
        )

    raise KeyError(
        f"""Baysor executable {default_path} does not exist. Please set a 'baysor' alias. Also check that you have installed baysor executable (as in https://github.com/kharchenkolab/Baysor)."""
    )
