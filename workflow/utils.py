from pathlib import Path
from typing import Optional


class WorkflowPaths:
    """
    A class listing the paths to different files that are needed
    or will be created by Snakemake (input and output)
    """

    def __init__(self, config: dict) -> None:
        self.config = sanity_check_config(config)

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

    def cells_paths(self, file_content: str, method_name: str):
        """Compute the paths to the temporary boundary files

        Args:
            file_content: Content of the file listing the number of patches or the patches indices
            method_name: Name of the method (cellpose, baysor, or comseg)
            dirs: Whether to return baysor directories

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

        # whether to run segmentation
        self.segmentation = "segmentation" in self.config  # whether to run segmentation

        # which segmentation method(s) is/are used
        self.cellpose = self.segmentation and "cellpose" in self.config["segmentation"]
        self.baysor = self.segmentation and "baysor" in self.config["segmentation"]
        self.comseg = self.segmentation and "comseg" in self.config["segmentation"]
        self.transcript_based_method = "comseg" if self.comseg else ("baysor" if self.baysor else None)

        # whether to run annotation
        self.annotate = "annotation" in self.config and "method" in self.config["annotation"]

        if self.baysor:
            check_baysor_executable_path(config)

    def resolve_transcripts(self) -> str:
        if self.transcript_based_method is None:
            return ""

        if "baysor" in self.config["segmentation"]:
            gene_column = self.config["segmentation"]["baysor"]["config"]["data"]["gene"]
        elif "comseg" in self.config["segmentation"]:
            gene_column = self.config["segmentation"]["comseg"]["config"]["gene_column"]

        min_area = self.config["segmentation"].get(self.transcript_based_method, {}).get("min_area", 0)
        return f"--gene-column {gene_column} --min-area {min_area}"

    def patchify_transcripts(self) -> str:
        if self.transcript_based_method is None:
            return ""

        dump_patch_params = str(self["patchify"].where(contains="micron"))
        method_config = self["segmentation"][self.transcript_based_method]

        if self.transcript_based_method == "comseg":
            dump_patch_params += " --write-cells-centroids"

        return f'{dump_patch_params} {method_config.where(keys=["prior_shapes_key", "unassigned_value"])}'

    ### The methods below are used to convert the Args object into a string for the Sopa CLI

    def __str__(self) -> str:
        """Convert an Args object into a string for the Sopa CLI"""
        return self.dump()

    def dump(self, prefix=""):
        return " ".join((res for item in self.config.items() for res in self.dump_arg(*item, prefix)))

    @classmethod
    def dump_arg(cls, key: str, value, prefix: str = ""):
        """
        Convert a key-value pair of the config into a string that
        can be used by the Sopa CLI
        """
        option = f"--{prefix}{key.replace('_', '-')}"
        if isinstance(value, list):
            for v in value:
                yield from (option, stringify_for_cli(v))
        elif value is True:
            yield option
        elif value is False:
            yield f"--no-{prefix}{key.replace('_', '-')}"
        else:
            yield from (option, stringify_for_cli(value))

    def __getitem__(self, name):
        subconfig = self.config.get(name, {})
        if not isinstance(subconfig, dict):
            return subconfig
        return Args(self.paths, subconfig)

    def where(self, keys: Optional[list[str]] = None, contains: Optional[str] = None):
        if keys is not None:
            return Args(self.paths, {key: self.config[key] for key in keys if key in self.config})
        if contains is not None:
            return Args(
                self.paths,
                {key: value for key, value in self.config.items() if contains in key},
            )
        return self


def stringify_for_cli(value) -> str:
    if isinstance(value, str) or isinstance(value, dict):
        return f'"{value}"'
    return str(value)


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


def sanity_check_config(config: dict):
    assert "segmentation" in config, "The `segmentation` parameter is mandatory"

    assert ("baysor" not in config["segmentation"]) or (
        "comseg" not in config["segmentation"]
    ), "Baysor and ComSeg cannot be used together"

    assert (
        "read" in config and "technology" in config["read"]
    ), "The `config['read']['technology'] parameter is mandatory"

    if config["read"]["technology"] in ["uniform", "toy_dataset"]:
        config["data_path"] = "."

    assert (
        "data_path" in config or "sdata_path" in config
    ), "Invalid config. Provide '--config data_path=...' when running the pipeline"

    if "data_path" in config and "sdata_path" not in config:
        config["sdata_path"] = Path(config["data_path"]).with_suffix(".zarr")
        print(
            f"SpatialData object path set to default: {config['sdata_path']}\nTo change this behavior, provide `--config sdata_path=...` when running the snakemake pipeline"
        )

    if "data_path" not in config:
        assert Path(
            config["sdata_path"]
        ).exists(), f"When `data_path` is not provided, the spatial data object must exist, but the directory doesn't exists: {config['sdata_path']}"
        config["data_path"] = []

    for method in ["baysor", "comseg"]:
        if method in config["segmentation"]:
            if "cellpose" in config["segmentation"]:
                config["segmentation"][method]["prior_shapes_key"] = "cellpose_boundaries"
            elif "cell_key" in config["segmentation"][method]:
                # backward compatibility
                config["segmentation"][method]["prior_shapes_key"] = config["segmentation"][method]["cell_key"]

    return config
