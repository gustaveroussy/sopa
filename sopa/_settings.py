from __future__ import annotations

import logging
import os
from typing import Callable

import dask
import dask.delayed
from dask.distributed import Client, progress
from tqdm import tqdm

dask.config.set({"dataframe.query-planning": False})  # SpatialData issue with dask-expr

log = logging.getLogger(__name__)


class Settings:
    ### General
    auto_save_on_disk: bool = True

    ### Parallelization
    _parallelization_backend: str | None = None
    available_parallelization_backends = [None, "dask"]
    dask_client_kwargs: dict = {}

    def __init__(self):
        self.parallelization_backend = os.environ.get("SOPA_PARALLELIZATION_BACKEND", None)

    @property
    def parallelization_backend(self):
        return self._parallelization_backend

    @parallelization_backend.setter
    def parallelization_backend(self, value):
        assert (
            value in self.available_parallelization_backends
        ), f"Invalid parallelization backend. Available options are: {self.available_parallelization_backends}"
        self._parallelization_backend = value

    def _run_with_backend(self, functions: list[Callable]):
        if self.parallelization_backend is None:
            if len(functions) > 1:
                log.warning(
                    "Running without parallelization backend can be slow. "
                    "Consider using a backend, e.g. via `sopa.settings.parallelization_backend = 'dask'`, "
                    "or `export SOPA_PARALLELIZATION_BACKEND=dask`."
                )
            return [f() for f in tqdm(functions)]

        log.info(f"Using {self.parallelization_backend} backend")
        return getattr(self, f"_run_{self.parallelization_backend}_backend")(functions)

    ### Dask backend
    def _run_dask_backend(self, functions: list[Callable]):
        assert len(functions) > 0, "No function to run"

        if len(functions) == 1:
            log.info("Only one function to run. No parallelization will be used.")
            functions[0]()

        import os

        @dask.delayed
        def run(f):
            return f()

        n_cpus = os.cpu_count()
        dask_client_kwargs = self.dask_client_kwargs.copy()
        if "n_workers" not in dask_client_kwargs and len(functions) < n_cpus:
            dask_client_kwargs["n_workers"] = n_cpus

        with Client(**dask_client_kwargs) as client:
            _ = dask.persist(*[run(f) for f in functions])
            progress(_, notebook=False)
            return client.gather(client.compute(_))


settings = Settings()
