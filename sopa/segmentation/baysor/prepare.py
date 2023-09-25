import toml


def to_toml(path: str, config: dict):
    with open(path, "w") as f:
        toml.dump(config, f)
