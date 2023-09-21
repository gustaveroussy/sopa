from pathlib import Path


class WorkflowPaths:
    def __init__(self, sdata_path: str) -> None:
        self.sdata_dir = Path(sdata_path)
        self.sdata = self.sdata_dir / ".zgroup"  # trick to fix snakemake ChildIOException
        self.raw = self.sdata_dir.with_suffix(".qptiff")  # TODO: make it general

        self.polygons = self.sdata_dir / "shapes" / "polygons"

        self.temp_dir = self.sdata_dir.parent / f"{self.sdata_dir.name}_temp"

        self.explorer_directory = self.sdata_dir.with_suffix(".explorer")
        self.explorer_directory.mkdir(parents=True, exist_ok=True)

        self.explorer_experiment = self.explorer_directory / "experiment.xenium"


def _dump_arg(key: str, value):
    option = f"--{key.replace('_', '-')}"
    if isinstance(value, list):
        for v in value:
            yield from (option, str(v))
    else:
        yield from (option, str(value))


def dump_args(args: dict) -> str:
    return " ".join((res for item in args.items() for res in _dump_arg(*item)))
