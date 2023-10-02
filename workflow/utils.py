from pathlib import Path


class ConfigConstants:
    CT_KEY = "cell_type"

    ANNOTATION = "annotation"


def _sanity_check_config(config: dict):
    for key in ["sdata_path"]:
        assert key in config.keys(), f"config['{key}'] is required to run the pipeline"

    return config


class WorkflowPaths:
    def __init__(self, config: dict) -> None:
        self.config = _sanity_check_config(config)

        self.sdata_path = Path(self.config["sdata_path"])
        self.sdata_zgroup = self.sdata_path / ".zgroup"  # trick to fix snakemake ChildIOException
        self.raw = self.sdata_path.with_suffix(".qptiff")  # TODO: make it general

        self.shapes_dir = self.sdata_path / "shapes"
        self.points_dir = self.sdata_path / "points"
        self.images_dir = self.sdata_path / "images"
        self.table_dir = self.sdata_path / "table"

        self.baysor_boundaries = self.shapes_dir / "baysor_boundaries"
        self.cellpose_boundaries = self.shapes_dir / "cellpose_boundaries"
        self.patches = self.shapes_dir / "patches"

        self.smk_files = self.sdata_path / ".smk_files"
        self.smk_table = self.smk_files / "table"
        self.smk_n_patches = self.smk_files / "n_patches"
        self.smk_aggregation = self.smk_files / "aggregation"

        self.annotations = []
        if "annotation" in self.config:
            self.annotations = (
                self.table_dir / "table" / "obs" / self.config["annotation"][ConfigConstants.CT_KEY]
            )

        self.temp_dir = self.sdata_path.parent / f"{self.sdata_path.name}_temp"
        self.cellpose_dir = self.temp_dir / "cellpose"
        self.baysor_dir = self.temp_dir / "baysor"

        self.explorer_directory = self.sdata_path.with_suffix(".explorer")
        self.explorer_directory.mkdir(parents=True, exist_ok=True)

        self.explorer_experiment = self.explorer_directory / "experiment.xenium"

    def cells_paths(self, n: int, name):
        if name == "cellpose":
            return [str(self.cellpose_dir / f"{i}.zarr.zip") for i in range(n)]
        if name == "baysor":
            return [str(self.baysor_dir / str(i) / "segmentation_polygons.json") for i in range(n)]


class Args:
    def __init__(self, paths: WorkflowPaths, config: dict):
        self.paths = paths
        self.config = config

        self.segmentation = "segmentation" in self.config
        self.cellpose = self.segmentation and "cellpose" in self.config["segmentation"]
        self.baysor = self.segmentation and "baysor" in self.config["segmentation"]
        self.annotate = ConfigConstants.ANNOTATION in self.config
        self.aggregate = "aggregate" in self.config

    def __getitem__(self, name):
        return Args(self.paths, self.config.get(name, {}))

    def dump_baysor_patchify(self):
        if not self.baysor:
            return ""

        return (
            self["segmentation"]["baysor"]._dump("baysor-")
            + f" --baysor-dir {self.paths.baysor_dir}"
        )

    @classmethod
    def _dump_arg(cls, key: str, value, prefix: str = ""):
        option = f"--{prefix}{key.replace('_', '-')}"
        if isinstance(value, list):
            for v in value:
                yield from (option, str(v))
        elif isinstance(value, dict):
            yield from (option, '"' + str(value) + '"')
        elif value is True:
            yield option
        elif value is False:
            yield f"--no-{prefix}{key.replace('_', '-')}"
        else:
            yield from (option, str(value))

    def _dump(self, prefix=""):
        return " ".join(
            (res for item in self.config.items() for res in self._dump_arg(*item, prefix))
        )

    def __str__(self) -> str:
        return self._dump()

    @property
    def baysor_cell_key(self):
        if not self.baysor or "cell_key" not in self.config["segmentation"]["baysor"]:
            return ""
        return f":{self.config['segmentation']['baysor']['cell_key']}"

    @property
    def gene_column(self):
        return self.config["segmentation"]["baysor"]["config"]["data"]["gene"]

    @property
    def expand_radius(self):
        return self.config["segmentation"]["cellpose"]["expand_radius"]
