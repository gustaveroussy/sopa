def _check_zip(names: list[str]):
    for name in names:
        if isinstance(name, str):
            assert name.endswith(".zip"), f"Intermediate files names must end with .zip"


SDATA_HELPER = "Path to the SpatialData `.zarr` directory"
