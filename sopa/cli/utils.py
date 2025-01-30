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


def _log_whether_to_resolve(patch_index: int | None, delete_cache: bool = True):
    import logging

    log = logging.getLogger(__name__)

    if patch_index is None and delete_cache:
        log.info("Segmentation is already resolved. Don't run `sopa resolve`.")


SDATA_HELPER = "Path to the SpatialData `.zarr` directory"
