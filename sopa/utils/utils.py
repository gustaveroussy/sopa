from __future__ import annotations

from pathlib import Path

from spatialdata import SpatialData

from .._constants import SopaFiles

HOME_CACHE_DIR = Path.home() / SopaFiles.SOPA_CACHE_DIR


def get_cache_dir(sdata: SpatialData) -> Path:
    cache_dir = sdata.path / SopaFiles.SOPA_CACHE_DIR if sdata.is_backed() else HOME_CACHE_DIR / str(id(sdata))

    cache_dir.mkdir(exist_ok=True, parents=True)

    return cache_dir


def clear_cache():
    import shutil

    for sub_dir in list(HOME_CACHE_DIR.iterdir()):
        shutil.rmtree(sub_dir)
