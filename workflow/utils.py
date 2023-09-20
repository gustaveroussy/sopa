from pathlib import Path


class WorkflowPaths:
    def __init__(self, sdata_path: str) -> None:
        self.sdata = Path(sdata_path)

        self.polygons = self.sdata / "shapes" / "polygons"

        self.temp_dir = self.sdata.parent / f".temp_{self.sdata.name}"


def _dump_arg(key: str, value):
    assert not isinstance(
        value, dict
    ), f"Can't dump dictionary {value} for key {key}. Instead, provide a str, int, float, or list"
    option = f"--{key.replace('_', '-')}"
    if isinstance(value, list):
        for v in value:
            yield from (option, v)
    else:
        yield from (option, value)


def dump_args(args: dict):
    return " ".join(map(str, (res for item in args.items() for res in _dump_arg(*item))))
