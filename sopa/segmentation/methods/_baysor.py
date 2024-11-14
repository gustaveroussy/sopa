import logging
from functools import partial
from pathlib import Path

from spatialdata import SpatialData

from ... import settings
from ..._constants import SopaAttrs, SopaFiles, SopaKeys
from ...utils import get_feature_key, get_transcripts_patches_dirs
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
    _check_transcript_patches(sdata)

    import shutil

    baysor_executable_path = _get_baysor_executable_path()
    use_polygons_format_argument = _use_polygons_format_argument(baysor_executable_path)

    prior_shapes_key = None
    if SopaKeys.PRIOR_SHAPES_KEY in sdata.shapes[SopaKeys.TRANSCRIPT_PATCHES]:
        prior_shapes_key = sdata.shapes[SopaKeys.TRANSCRIPT_PATCHES][SopaKeys.PRIOR_SHAPES_KEY].iloc[0]

    if config is None or not len(config):
        config = _get_default_config(sdata, prior_shapes_key, scale)

    patches_dirs = get_transcripts_patches_dirs(sdata)
    for patch_dir in patches_dirs:
        _copy_segmentation_config(patch_dir / SopaFiles.TOML_CONFIG_FILE, config)

    gene_column = _get_gene_column_argument(config)

    baysor_patch = BaysorPatch(
        baysor_executable_path,
        use_polygons_format_argument,
        force=force,
        recover=recover,
        prior_shapes_key=prior_shapes_key,
    )

    if patch_index is not None:
        baysor_patch(patches_dirs[patch_index])
        return

    settings._run_with_backend([partial(baysor_patch, patch_dir) for patch_dir in patches_dirs])

    if force:
        patches_dirs = [patch_dir for patch_dir in patches_dirs if (patch_dir / "segmentation_counts.loom").exists()]
        assert patches_dirs, "Baysor failed on all patches"

    resolve(sdata, patches_dirs, gene_column, min_area=min_area, key_added=key_added)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added

    if delete_cache:
        for patch_dir in get_transcripts_patches_dirs(sdata):
            shutil.rmtree(patch_dir)


class BaysorPatch:
    def __init__(
        self,
        baysor_executable_path: str,
        use_polygons_format_argument: bool,
        force: bool = False,
        recover: bool = False,
        prior_shapes_key: str | None = None,
    ):
        self.baysor_executable_path = baysor_executable_path
        self.use_polygons_format_argument = use_polygons_format_argument
        self.force = force
        self.recover = recover
        self.prior_shapes_key = prior_shapes_key

    def __call__(self, patch_dir: Path):
        if self.recover and (patch_dir / "segmentation_counts.loom").exists():
            return

        import subprocess

        polygon_substring = (
            "--polygon-format GeometryCollection" if self.use_polygons_format_argument else "--save-polygons GeoJSON"
        )

        prior_suffix = f":{self.prior_shapes_key}" if self.prior_shapes_key else ""

        baysor_command = (
            f"{self.baysor_executable_path} run {polygon_substring} -c config.toml transcripts.csv {prior_suffix}"
        )

        result = subprocess.run(
            f"""
            cd {patch_dir}
            {baysor_command}
            """,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if result.returncode != 0:
            message = f"Baysor error on patch {patch_dir} with command `{baysor_command}`"
            if self.force:
                log.warning(message)
                return
            raise RuntimeError(f"{message}:\n{result.stdout.decode()}")


def _use_polygons_format_argument(baysor_executable_path: str) -> bool:
    import subprocess

    from packaging.version import Version

    try:
        res = subprocess.run(
            f"{baysor_executable_path} --version",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        return Version(res.stdout) >= Version("0.7.0")
    except:
        return False


def _get_baysor_executable_path() -> Path | str:
    import shutil

    if shutil.which("baysor") is not None:
        return "baysor"

    default_path = Path.home() / ".julia" / "bin" / "baysor"
    if default_path.exists():
        return default_path

    raise FileNotFoundError(
        f"Please install baysor and ensure that either `{default_path}` executes baysor, or `baysor` is an existing shell alias for baysor's executable."
    )


def _get_default_config(sdata: SpatialData, prior_shapes_key: str | None, scale: float | None) -> dict:
    assert prior_shapes_key is not None or scale is not None, (
        "The config can't be inferred. Choose among the following solutions:\n"
        "   - Provide the `scale` argument (typical cell radius in microns)\n"
        "   - Run `sopa.make_transcript_patches(...)` with `prior_shapes_key`\n"
        "   - Provide the `config` argument, containing a valid Baysor config."
    )

    points_key = sdata.attrs.get(SopaAttrs.TRANSCRIPTS)
    assert (
        points_key
    ), f"Transcripts key not found in sdata.attrs['{SopaAttrs.TRANSCRIPTS}'], baysor config can't be inferred."

    feature_key = get_feature_key(sdata[points_key], raise_error=True)

    config = {
        "data": {
            "x": "x",
            "y": "y",
            "gene": feature_key,
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
