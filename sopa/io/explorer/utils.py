from pathlib import Path


def explorer_file_path(path: str, filename: str, is_dir: bool):
    path: Path = Path(path)

    if not is_dir:
        path = path / filename

    return path
