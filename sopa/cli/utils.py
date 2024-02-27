from __future__ import annotations


def _check_zip(names: list[str]):
    for name in names:
        if isinstance(name, str):
            assert name.endswith(".zip"), "Intermediate files names must end with .zip"


def _default_boundary_dir(sdata_path: str, directory_name: str):
    from pathlib import Path

    from sopa._constants import SopaFiles

    temp_dir = Path(sdata_path) / SopaFiles.SOPA_CACHE_DIR / directory_name
    temp_dir.mkdir(parents=True, exist_ok=True)

    return temp_dir


SDATA_HELPER = "Path to the SpatialData `.zarr` directory"
