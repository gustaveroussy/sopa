try:
    import toml
except ImportError:
    raise ImportError(
        "To use baysor, you need its corresponding sopa extra: `pip install 'sopa[baysor]'`.\
        \nAlso, make sure to install the baysor executable (https://github.com/kharchenkolab/Baysor)."
    )


def to_toml(path: str, config: dict):
    with open(path, "w") as f:
        toml.dump(config, f)
