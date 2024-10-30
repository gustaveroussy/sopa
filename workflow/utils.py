from __future__ import annotations

from pathlib import Path
from typing import Optional


class WorkflowPaths:
    """
    A class listing the paths to different files that are needed
    or will be created by Snakemake (input and output)
    """

    def __init__(self, config: dict) -> None:
        self.config = config

        ### SpatialData object files
        self.data_path = self.config["data_path"]
        self.sdata_path = Path(self.config["sdata_path"])
        self.sdata_zgroup = self.sdata_path / ".zgroup"  # trick to fix snakemake ChildIOException

        self.shapes_dir = self.sdata_path / "shapes"
        self.points_dir = self.sdata_path / "points"
        self.images_dir = self.sdata_path / "images"
        self.table_dir = self.sdata_path / "tables"

        self.baysor_boundaries = self.shapes_dir / "baysor_boundaries"
        self.cellpose_boundaries = self.shapes_dir / "cellpose_boundaries"

        ### Snakemake cache / intermediate files
        self.sopa_cache = self.sdata_path / ".sopa_cache"

        self.smk_patches = self.sopa_cache / "patches"
        self.smk_aggregation = self.sopa_cache / "aggregation"
        self.smk_table = self.sopa_cache / "table"

        self.smk_cellpose_temp_dir = self.sopa_cache / "cellpose_boundaries"
        self.smk_patches_file_image = self.sopa_cache / "patches_file_image"

        self.smk_transcripts_temp_dir = self.sopa_cache / "transcript_patches"
        self.smk_patches_file_transcripts = self.sopa_cache / "patches_file_transcripts"

        self.smk_cellpose_boundaries = self.sopa_cache / "cellpose_boundaries_done"
        self.smk_baysor_boundaries = self.sopa_cache / "baysor_boundaries_done"
        self.smk_comseg_boundaries = self.sopa_cache / "comseg_boundaries_done"

        ### Xenium Explorer outputs
        self.explorer_directory = self.sdata_path.with_suffix(".explorer")
        self.explorer_directory.mkdir(parents=True, exist_ok=True)
        self.explorer_experiment = self.explorer_directory / "experiment.xenium"
        self.explorer_image = self.explorer_directory / "morphology.ome.tif"
        self.report = self.explorer_directory / "analysis_summary.html"

        ### Annotation files
        self.annotations = []
        if "annotation" in self.config:
            key = self.config["annotation"].get("args", {}).get("cell_type_key", "cell_type")
            self.annotations = self.table_dir / "table" / "obs" / key

    def temporary_boundaries_paths(self, file_content: str, method_name: str) -> list[str]:
        """Compute the paths to the temporary boundary files

        Args:
            file_content: Content of the file listing the number of patches (or the patches indices)
            method_name: Name of the method: cellpose, baysor, or comseg

        Returns:
            A list of temporary boundary directories or files
        """
        if method_name == "cellpose":
            return [str(self.smk_cellpose_temp_dir / f"{i}.parquet") for i in range(int(file_content))]

        if method_name == "baysor":
            indices = map(int, file_content.split())
            return [str(self.smk_transcripts_temp_dir / str(i) / "segmentation_counts.loom") for i in indices]

        if method_name == "comseg":
            indices = map(int, file_content.split())
            COMSEG_FILES = ["segmentation_polygons.json", "segmentation_counts.h5ad"]
            return [str(self.smk_transcripts_temp_dir / str(i) / file) for i in indices for file in COMSEG_FILES]


class Args:
    """
    A convenient class to pass the YAML config arguments to sopa's CLI
    """

    def __init__(self, paths: WorkflowPaths, config: dict):
        self.paths = paths
        self.config = config

        # which segmentation method(s) is/are used
        self.cellpose = "cellpose" in self.config.get("segmentation", {})
        self.baysor = "baysor" in self.config.get("segmentation", {})
        self.comseg = "comseg" in self.config.get("segmentation", {})
        self.transcript_based_method = "comseg" if self.comseg else ("baysor" if self.baysor else None)

        # whether to run annotation
        self.annotate = "annotation" in self.config and "method" in self.config["annotation"]

    def resolve_transcripts(self) -> str:
        """Arguments for `sopa resolve [baysor/comseg]`"""
        if self.transcript_based_method is None:
            return ""

        if "baysor" in self.config["segmentation"]:
            gene_column = self.config["segmentation"]["baysor"]["config"]["data"]["gene"]
        elif "comseg" in self.config["segmentation"]:
            gene_column = self.config["segmentation"]["comseg"]["config"]["gene_column"]

        min_area = self.config["segmentation"].get(self.transcript_based_method, {}).get("min_area", 0)
        return f"--gene-column {gene_column} --min-area {min_area}"

    def patchify_transcripts(self) -> str:
        """Arguments for `sopa patchify transcripts`"""
        if self.transcript_based_method is None:
            return ""

        params = self["patchify"].as_cli(contains="micron")

        if self.transcript_based_method == "comseg":
            params += " --write-cells-centroids"

        method_config = self["segmentation"][self.transcript_based_method]
        return f'{params} {method_config.as_cli(keys=["prior_shapes_key", "unassigned_value"])}'

    ### The methods below are used to convert the Args object into a string for the Sopa CLI

    def as_cli(self, keys: Optional[list[str]] = None, contains: Optional[str] = None) -> str:
        """Extract a subset of the config (or the whole config) as a string for the CLI (command-line interface)

        Args:
            keys: List of keys to extract from the config.
            contains: String that must be contained in the keys to be extracted.

        Returns:
            A string that can be used as arguments/options for the Sopa CLI.
        """
        assert (keys is None) or (contains is None), "Provide either 'keys' or 'contains', but not both"

        if keys is None and contains is None:
            return str(self)

        if keys is not None:
            sub_args = Args(self.paths, {key: self.config[key] for key in keys if key in self.config})
        elif contains is not None:
            sub_args = Args(self.paths, {key: value for key, value in self.config.items() if contains in key})

        return str(sub_args)

    def __str__(self) -> str:
        """
        Config as a string, useful for the Sopa CLI

        For instance, {"x": 2, "y": False} will be converted to "--x 2 --no-y"
        """
        return " ".join((res for item in self.config.items() for res in stringify_item(*item)))

    def __getitem__(self, name: str) -> Args | bool | str | list:
        sub_config = self.config.get(name, {})
        if not isinstance(sub_config, dict):
            return sub_config
        return Args(self.paths, sub_config)


def stringify_item(key: str, value: bool | list | dict | str):
    """
    Convert a key-value pair of the config into a string that can be used by the Sopa CLI
    """
    key = key.replace("_", "-")
    option = f"--{key}"

    if value is True:
        yield option
    elif value is False:
        yield f"--no-{key}"
    elif isinstance(value, list):
        for v in value:
            yield from (option, _stringify_value_for_cli(v))
    else:
        yield from (option, _stringify_value_for_cli(value))


def _stringify_value_for_cli(value) -> str:
    if isinstance(value, str) or isinstance(value, dict):
        return f'"{value}"'
    return str(value)
