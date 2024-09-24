# Adapted from https://github.com/manzt/napari-lazy-openslide/tree/main
from ctypes import ArgumentError
from pathlib import Path
from typing import Any, Dict, Mapping, MutableMapping

import numpy as np
from openslide import OpenSlide
from zarr.storage import (
    KVStore,
    Store,
    _path_to_prefix,
    attrs_key,
    init_array,
    init_group,
)
from zarr.util import json_dumps, normalize_shape, normalize_storage_path


def init_attrs(store: MutableMapping, attrs: Mapping[str, Any], path: str = None):
    path = normalize_storage_path(path)
    path = _path_to_prefix(path)
    store[path + attrs_key] = json_dumps(attrs)


def create_meta_store(slide: OpenSlide, tilesize: int) -> Dict[str, bytes]:
    """Creates a dict containing the zarr metadata for the multiscale openslide image."""
    store = dict()
    root_attrs = {
        "multiscales": [
            {
                "name": Path(slide._filename).name,
                "datasets": [{"path": str(i)} for i in range(slide.level_count)],
                "version": "0.1",
            }
        ],
        "metadata": dict(slide.properties),
    }
    init_group(store)
    init_attrs(store, root_attrs)

    for i, (x, y) in enumerate(slide.level_dimensions):
        init_array(
            store,
            path=str(i),
            shape=normalize_shape((y, x, 4)),
            chunks=(tilesize, tilesize, 4),
            fill_value=0,
            dtype="|u1",
            compressor=None,
        )
        suffix = str(i) if i != 0 else ""
        init_attrs(store, {"_ARRAY_DIMENSIONS": [f"Y{suffix}", f"X{suffix}", "S"]}, path=str(i))
    return store


def _parse_chunk_path(path: str):
    """Returns x,y chunk coords and pyramid level from string key"""
    level, ckey = path.split("/")
    y, x, _ = map(int, ckey.split("."))
    return x, y, int(level)


class OpenSlideStore(Store):
    """Wraps an OpenSlide object as a multiscale Zarr Store.

    Parameters
    ----------
    path: str
        The file to open with OpenSlide.
    tilesize: int
        Desired "chunk" size for zarr store (default: 512).
    """

    def __init__(self, path: str, tilesize: int = 512):
        self._path = path
        self._slide = OpenSlide(path)
        self._tilesize = tilesize
        self._store = create_meta_store(self._slide, tilesize)
        self._writeable = False
        self._erasable = False

    def __getitem__(self, key: str):
        if key in self._store:
            # key is for metadata
            return self._store[key]

        # key should now be a path to an array chunk
        # e.g '3/4.5.0' -> '<level>/<chunk_key>'
        try:
            x, y, level = _parse_chunk_path(key)
            location = self._ref_pos(x, y, level)
            size = (self._tilesize, self._tilesize)
            tile = self._slide.read_region(location, level, size)
        except ArgumentError as err:
            # Can occur if trying to read a closed slide
            raise err
        except Exception:
            # TODO: probably need better error handling.
            # If anything goes wrong, we just signal the chunk
            # is missing from the store.
            raise KeyError(key)
        return np.array(tile)  # .tobytes()

    def __eq__(self, other):
        return isinstance(other, OpenSlideStore) and self._slide._filename == other._slide._filename

    def __setitem__(self, key, val):
        raise PermissionError("ZarrStore is read-only")

    def __delitem__(self, key):
        raise PermissionError("ZarrStore is read-only")

    def __iter__(self):
        return iter(self.keys())

    def __len__(self):
        return sum(1 for _ in self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _ref_pos(self, x: int, y: int, level: int):
        dsample = self._slide.level_downsamples[level]
        xref = int(x * dsample * self._tilesize)
        yref = int(y * dsample * self._tilesize)
        return xref, yref

    def keys(self):
        return self._store.keys()

    def close(self):
        self._slide.close()

    # Retrieved from napari-lazy-openslide PR#16
    def __getstate__(self):
        return {"_path": self._path, "_tilesize": self._tilesize}

    def __setstate__(self, newstate):
        path = newstate["_path"]
        tilesize = newstate["_tilesize"]
        self.__init__(path, tilesize)

    def rename(self):
        raise PermissionError(f'{type(self)} is not erasable, cannot call "rename"')

    def rmdir(self, path: str = "") -> None:
        raise PermissionError(f'{type(self)} is not erasable, cannot call "rmdir"')

    @property
    def store(self) -> KVStore:
        return KVStore(self)
