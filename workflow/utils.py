from pathlib import Path


class WorkflowPaths:
    def __init__(self, sdata_path: str, config_path: str) -> None:
        self.sdata_path = Path(sdata_path)
        self.config_path = Path(config_path)

        self.temp_dir = self.sdata_path.parent / f".temp_{self.sdata_path.name}"
