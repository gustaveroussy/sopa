import shutil
from pathlib import Path


def _get_executable_path(name: str, default_dir: str) -> Path | str:
    if shutil.which(name) is not None:
        return name

    default_path = Path.home() / default_dir / "bin" / name
    if default_path.exists():
        return default_path

    bin_path = Path.home() / ".local" / "bin" / name
    raise FileNotFoundError(
        f"Please install {name} and ensure that either `{default_path}` executes {name}, or that `{name}` is an existing command (add it to your PATH, or create a symlink at {bin_path})."
    )
