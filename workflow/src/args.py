from __future__ import annotations

from .constants import SEGMENTATION_METHODS, TRANSCRIPT_BASED_METHODS
from .paths import WorkflowPaths


class Args:
    """
    A convenient class to pass the YAML config arguments to sopa's CLI
    """

    def __init__(self, paths: WorkflowPaths, config: dict):
        self.paths = paths
        self.config = config

        # which transcript-based segmentation to run (if any)
        self.transcript_based_method = None
        for method in TRANSCRIPT_BASED_METHODS:
            if method in self.config.get("segmentation", {}):
                self.transcript_based_method = method
                break

        # whether to run annotation
        self.annotate = "annotation" in self.config and "method" in self.config["annotation"]

    def use(self, method_name: str) -> bool:
        return method_name in self.config["segmentation"]

    def segmentation_boundaries(self):
        for method in SEGMENTATION_METHODS:
            if self.use(method):
                return self.paths.segmentation_done(method)
        raise ValueError("No segmentation method selected")

    def resolve_transcripts(self) -> str:
        """Arguments for `sopa resolve [baysor/comseg]`"""
        if self.transcript_based_method is None or self.transcript_based_method == "proseg":
            return ""

        if self.transcript_based_method == "baysor":
            gene_column = self.config["segmentation"]["baysor"]["config"]["data"]["gene"]
        elif self.transcript_based_method == "comseg":
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
        return f"{params} {method_config.as_cli(keys=['prior_shapes_key', 'unassigned_value'])}"

    ### The methods below are used to convert the Args object into a string for the Sopa CLI

    def as_cli(self, keys: list[str] | None = None, contains: str | None = None) -> str:
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
            sub_args = Args(
                self.paths,
                {key: self.config[key] for key in keys if key in self.config},
            )
        elif contains is not None:
            sub_args = Args(
                self.paths,
                {key: value for key, value in self.config.items() if contains in key},
            )

        return str(sub_args)

    def __str__(self) -> str:
        """
        Config as a string, useful for the Sopa CLI

        For instance, {"x": 2, "y": False} will be converted to "--x 2 --no-y"
        """
        return " ".join(res for item in self.config.items() for res in _stringify_item(*item))

    def __getitem__(self, name: str) -> Args | bool | str | list:
        sub_config = self.config.get(name, {})
        if not isinstance(sub_config, dict):
            return sub_config
        return Args(self.paths, sub_config)


def _stringify_item(key: str, value: bool | list | dict | str):
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
    if isinstance(value, (str, dict)):
        return f'"{value}"'
    return str(value)
