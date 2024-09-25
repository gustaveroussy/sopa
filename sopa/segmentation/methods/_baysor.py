from __future__ import annotations

import logging
from pathlib import Path

from spatialdata import SpatialData

from ..._constants import SopaFiles, SopaKeys
from ..transcripts import copy_segmentation_config, resolve

log = logging.getLogger(__name__)


def baysor(sdata: SpatialData, config: dict | None = None, config_path: str | None = None, min_area: int = 0):
    import shutil

    baysor_executable_path = _get_baysor_executable_path()
    use_polygons_format_argument = _use_polygons_format_argument(baysor_executable_path)

    assert (config is None) ^ (config_path is None), "Provide either a config dict or a path to a config file"

    if config_path is not None:
        import toml

        config = toml.load(config_path)

    assert config.get("data", {}).get("gene"), "Gene column not found in config['data']['gene']"
    gene_column = config["data"]["gene"]

    assert (
        SopaKeys.TRANSCRIPT_PATCHES in sdata.shapes
    ), "Transcript patches not found in the SpatialData object. Run `sopa.make_transcript_patches(...)` first."

    patches_dirs = sdata.shapes[SopaKeys.TRANSCRIPT_PATCHES][SopaKeys.CACHE_PATH_KEY].values

    for patch_dir in patches_dirs:
        copy_segmentation_config(Path(patch_dir) / SopaFiles.TOML_CONFIG_FILE, config, config_path)
        baysor_patch(patch_dir, baysor_executable_path, use_polygons_format_argument)

    resolve(sdata, None, gene_column, min_area=min_area, patches_dirs=patches_dirs)

    for patch_dir in patches_dirs:
        shutil.rmtree(patch_dir)


def baysor_patch(patch_dir: str, baysor_executable_path: str, use_polygons_format_argument: bool):
    import subprocess

    polygon_substring = (
        "--polygon-format GeometryCollection" if use_polygons_format_argument else "--save-polygons GeoJSON"
    )

    command = f"""
        cd {patch_dir}
        {baysor_executable_path} run {polygon_substring} -c config.toml transcripts.csv
        """
    subprocess.run(command, shell=True)


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
