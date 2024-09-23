import logging
from typing import Callable

import dask
import dask.delayed
from dask.distributed import Client
from tqdm import tqdm

dask.config.set({"dataframe.query-planning": False})  # SpatialData issue with dask-expr

log = logging.getLogger(__name__)


class Settings:
    auto_save_on_disk: bool = True
    AVAILABLE_PARALLELIZATION_BACKENDS = [None, "dask"]

    def __init__(self) -> None:
        self._parallelization_backend = None

    @property
    def parallelization_backend(self):
        return self._parallelization_backend

    @parallelization_backend.setter
    def parallelization_backend(self, value):
        assert (
            value in self.AVAILABLE_PARALLELIZATION_BACKENDS
        ), f"Invalid parallelization backend. Available options are: {self.AVAILABLE_PARALLELIZATION_BACKENDS}"
        self._parallelization_backend = value

    def _run_with_backend(self, functions: list[Callable]):
        if self.parallelization_backend is None:
            log.warning(
                "Running without parallelization backend can be slow. Consider using a backend, e.g. `sopa.settings.parallelization_backend = 'dask'`."
            )
            return [f() for f in tqdm(functions)]
        else:
            log.info(f"Using {self.parallelization_backend} backend")
            return getattr(self, f"_run_{self.parallelization_backend}_backend")(functions)

    ### Dask backend
    def _run_dask_backend(self, functions: list[Callable]):
        client = Client()

        @dask.delayed
        def run(f):
            return f()

        res = dask.compute(*[run(f) for f in functions])
        client.close()

        return res


settings = Settings()
