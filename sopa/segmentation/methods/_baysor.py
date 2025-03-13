import logging
from functools import partial
from pathlib import Path

from spatialdata import SpatialData

from ... import settings
from ..._constants import SopaAttrs, SopaFiles, SopaKeys
from ...utils import (
    delete_transcripts_patches_dirs,
    get_feature_key,
    get_transcripts_patches_dirs,
)
from .._transcripts import _check_transcript_patches, resolve

log = logging.getLogger(__name__)


def baysor(
    sdata: SpatialData,
    config: dict | str | None = None,
    min_area: int = 0,
    delete_cache: bool = True,
    recover: bool = False,
    force: bool = False,
    scale: float | None = None,
    key_added: str = SopaKeys.BAYSOR_BOUNDARIES,
    patch_index: int | None = None,
):
    """Run [Baysor](https://kharchenkolab.github.io/Baysor/dev/) segmentation on a SpatialData object, and add a GeoDataFrame containing the cell boundaries.

    !!! warning "Baysor installation"
        Make sure to install [Baysor](https://kharchenkolab.github.io/Baysor/dev/installation/), and either have the executable at `~/.julia/bin/baysor`, or create an alias called
        `baysor` that points to the binary executable. Also, you'll need to install
        sopa with the baysor extra: `pip install 'sopa[baysor]'` (basically, this installs `toml` and `loompy`).

    !!! info "Inferred config"
        If the `config` argument is not provided, the configuration is inferred.
        If [sopa.make_transcript_patches][] was run with a `prior_shapes_key`, the configuration is inferred based on the prior segmentation.
        Otherwise, the configuration is inferred based on the `scale` parameter (you'll need to provide it).

    Args:
        sdata: A `SpatialData` object.
        config: Optional configuration dictionary or path to a TOML file containing a valid Baysor config. By default, a configuration is inferred based on the cell area of the prior segmentation, or based on the `scale` parameter.
        min_area: Minimal area (in microns^2) of a cell to be considered.
        delete_cache: Whether to delete the cache after segmentation.
        recover: If `True`, recover the cache from a failed segmentation, and continue.
        force: If `True`, ignore failed patches and continue with the successful ones.
        scale: The typical cell radius in microns. If `config` is not provided, the configuration is inferred based on this parameter.
        key_added: Name of the shapes element to be added to `sdata.shapes`.
        patch_index: Index of the patch to segment (we do not recommend to set this argument). By default, segment all patches.
    """
    _check_transcript_patches(sdata)

    prior_shapes_key = None
    if SopaKeys.PRIOR_SHAPES_KEY in sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES]:
        prior_shapes_key = sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.PRIOR_SHAPES_KEY].iloc[0]

    if config is None or not len(config):
        config = _get_default_config(sdata, prior_shapes_key, scale)

    baysor_command = _get_baysor_command(prior_shapes_key)

    baysor_patch = BaysorPatch(baysor_command, config, force=force, capture_output=patch_index is None)

    if patch_index is not None:
        patch_dir = Path(sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES].loc[patch_index, SopaKeys.CACHE_PATH_KEY])
        baysor_patch(patch_dir)
        return

    patches_dirs = get_transcripts_patches_dirs(sdata)

    remaining_patches_dirs = (
        [patch_dir for patch_dir in patches_dirs if not (patch_dir / "segmentation_counts.loom").exists()]
        if recover
        else patches_dirs
    )

    settings._run_with_backend([partial(baysor_patch, patch_dir) for patch_dir in remaining_patches_dirs])

    if force:
        patches_dirs = [patch_dir for patch_dir in patches_dirs if (patch_dir / "segmentation_counts.loom").exists()]
        assert patches_dirs, "Baysor failed on all patches"

    gene_column = _get_gene_column_argument(config)
    resolve(sdata, patches_dirs, gene_column, min_area=min_area, key_added=key_added)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added

    if delete_cache:
        delete_transcripts_patches_dirs(sdata)


class BaysorPatch:
    def __init__(
        self,
        baysor_command: str,
        config: dict | str,
        force: bool = False,
        capture_output: bool = True,
    ):
        self.baysor_command = baysor_command
        self.config = config
        self.force = force
        self.capture_output = capture_output

    def __call__(self, patch_dir: Path):
        _copy_segmentation_config(patch_dir / SopaFiles.TOML_CONFIG_FILE, self.config)

        import subprocess

        result = subprocess.run(self.baysor_command, cwd=patch_dir, shell=True, capture_output=self.capture_output)

        if result.returncode != 0 or not (patch_dir / "segmentation_counts.loom").exists():
            if self.force:
                log.warning(f"Baysor error on patch {patch_dir.resolve()} with command `{self.baysor_command}`")
                return
            raise subprocess.CalledProcessError(
                returncode=result.returncode,
                cmd=self.baysor_command,
                output=result.stdout,
                stderr=result.stderr,
            )


def _get_baysor_command(prior_shapes_key: str | None) -> str:
    baysor_executable_path = _get_baysor_executable_path()

    use_polygons_format_argument = _use_polygons_format_argument(baysor_executable_path)
    polygon_format = (
        "--polygon-format GeometryCollection" if use_polygons_format_argument else "--save-polygons GeoJSON"
    )  # depends on the version of baysor

    prior_suffix = f" :{prior_shapes_key}" if prior_shapes_key else ""  # use a prior segmentation

    return f"{baysor_executable_path} run {polygon_format} -c config.toml transcripts.csv" + prior_suffix


def _use_polygons_format_argument(baysor_executable_path: str) -> bool:
    import subprocess

    from packaging.version import InvalidVersion, Version

    result = subprocess.run(
        f"{baysor_executable_path} run --version",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    try:
        return Version(result.stdout) >= Version("0.7.0")
    except InvalidVersion:
        log.warning("Could not parse the version of baysor. Assumes baysor >= 0.7.0.")
        return True


def _get_baysor_executable_path() -> Path | str:
    import shutil

    if shutil.which("baysor") is not None:
        return "baysor"

    default_path = Path.home() / ".julia" / "bin" / "baysor"
    if default_path.exists():
        return default_path

    bin_path = Path.home() / ".local" / "bin" / "baysor"
    raise FileNotFoundError(
        f"Please install baysor and ensure that either `{default_path}` executes baysor, or that `baysor` is an existing command (add it to your PATH, or create a symlink at {bin_path})."
    )


def _get_default_config(sdata: SpatialData, prior_shapes_key: str | None, scale: float | None) -> dict:
    assert prior_shapes_key is not None or scale is not None, (
        "The config can't be inferred. Choose among the following solutions:\n"
        "   - Provide the `scale` argument (typical cell radius in microns)\n"
        "   - Run `sopa.make_transcript_patches(...)` with `prior_shapes_key`\n"
        "   - Provide the `config` argument, containing a valid Baysor config."
    )

    points_key = sdata[SopaKeys.TRANSCRIPTS_PATCHES][SopaKeys.POINTS_KEY].iloc[0]
    feature_key = get_feature_key(sdata[points_key], raise_error=True)

    config = {
        "data": {
            "x": "x",
            "y": "y",
            "gene": str(feature_key),
            "min_molecules_per_gene": 10,
            "min_molecules_per_cell": 20,
            "force_2d": True,
        },
        "segmentation": {"prior_segmentation_confidence": 0.8},
    }

    if scale is not None:
        config["segmentation"]["scale"] = scale

    log.info(f"The Baysor config was not provided, using the following by default:\n{config}")

    return config


def _get_gene_column_argument(config: dict | str) -> str:
    if isinstance(config, str):
        import toml

        config = toml.load(config)

    assert config.get("data", {}).get("gene"), "Gene column not found in config['data']['gene']"
    return config["data"]["gene"]


def _copy_segmentation_config(path: Path | str, config: dict | str):
    """Copy the segmentation config to a file (`.json` or `.toml`).

    Args:
        path: Where the config will be saved
        config: Dictionnary config, or path to an existing config file (json or toml)
    """
    path = Path(path)

    if isinstance(config, str):
        import shutil

        shutil.copy(config, path)
        return

    assert path.suffix == ".toml"

    try:
        import toml
    except ImportError:
        raise ImportError(
            "To use baysor, you need its corresponding sopa extra: `pip install 'sopa[baysor]'`.\
            \nAlso, make sure to install the baysor executable (https://github.com/kharchenkolab/Baysor)."
        )

    with open(path, "w") as f:
        toml.dump(config, f)
        return
