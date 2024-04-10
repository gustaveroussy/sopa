from __future__ import annotations

try:
    import toml
except ImportError:
    raise ImportError(
        "To use baysor, you need its corresponding sopa extra: `pip install 'sopa[baysor]'` (normal mode) or `pip install -e '.[baysor]'` (if using snakemake).\
        \nAlso, make sure to install the baysor executable (https://github.com/kharchenkolab/Baysor)."
    )


def copy_toml_config(path: str, config: dict, config_path: str | None):
    if config_path is not None:
        import shutil

        shutil.copy(config_path, path)
    else:
        with open(path, "w") as f:
            toml.dump(config, f)
