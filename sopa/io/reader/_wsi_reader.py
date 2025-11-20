import logging

from xarray import DataArray, DataTree

log = logging.getLogger(__name__)


class ReaderBase:
    path: str
    name = "base"
    slide = None

    def read_region(self, location, level, size, **kwargs):
        raise NotImplementedError

    def get_metadata(self):
        """Get metadata from the slide."""
        return {
            "properties": self.properties,
            "dimensions": self.dimensions,
            "level_count": self.level_count,
            "level_dimensions": self.level_dimensions,
            "level_downsamples": self.level_downsamples,
        }

    def get_zarr_store(self):
        """Get the Zarr store for the slide."""
        raise NotImplementedError

    def close(self):
        raise NotImplementedError

    @property
    def properties(self):
        raise NotImplementedError

    @property
    def level_downsamples(self):
        raise NotImplementedError

    @property
    def level_count(self):
        raise NotImplementedError

    @property
    def level_dimensions(self):
        raise NotImplementedError

    @property
    def dimensions(self):
        raise NotImplementedError


class OpenSlideReader(ReaderBase):
    """Reader using the openslide backend."""

    name = "openslide"

    def __init__(self, path: str):
        try:
            from openslide import OpenSlide
        except ImportError:
            raise ImportError(
                "To use the openslide backend, you need to install it, e.g.,: `pip install openslide-python openslide-bin`."
            )

        self.path = path
        self.slide = OpenSlide(path)

    def read_region(self, location, level, size):
        """Get a region from the slide."""
        return self.slide.read_region(location, level, size).convert("RGB")

    def get_zarr_store(self, tilesize: int = 512):
        from zarr.storage import KVStore

        from ._wsi_store import WsiStore

        return KVStore(WsiStore(self, tilesize))

    def close(self):
        """Close the slide."""
        self.slide.close()

    @property
    def properties(self):
        """Get the properties of the slide."""
        return self.slide.properties

    @property
    def level_downsamples(self):
        return self.slide.level_downsamples

    @property
    def level_count(self):
        return len(self.slide.level_downsamples)

    @property
    def level_dimensions(self):
        return self.slide.level_dimensions

    @property
    def dimensions(self):
        """Get the dimensions of the slide."""
        return self.slide.dimensions


class TiffSlideReader(ReaderBase):
    """Reader using the tiffslide backend."""

    name = "tiffslide"

    def __init__(self, path: str):
        import tiffslide

        self.path = path
        self.slide = tiffslide.open_slide(path)

    def read_region(self, location, level, size):
        """Get a region from the slide."""
        return self.slide.read_region(location, level, size).convert("RGB")

    def get_zarr_store(self, tilesize: int = 512):
        return self.slide.zarr_group.store

    def close(self):
        """Close the slide."""
        self.slide.close()

    @property
    def properties(self):
        """Get the properties of the slide."""
        return self.slide.properties

    @property
    def level_downsamples(self):
        return self.slide.level_downsamples

    @property
    def level_count(self):
        return len(self.slide.level_downsamples)

    @property
    def level_dimensions(self):
        return self.slide.level_dimensions

    @property
    def dimensions(self):
        """Get the dimensions of the slide."""
        return self.slide.dimensions


class SlideIOReader(ReaderBase):
    """Reader using the SlideIO backend."""

    name = "slideio"

    def __init__(self, path: str):
        import dask

        dask.config.set(scheduler="single-threaded")
        log.warning("SlideIOReader is not multi-threaded compatible, setting dask scheduler to single-threaded.")

        try:
            import slideio
        except ImportError:
            raise ImportError("To use the slideio backend, you need to install it, e.g.,: `pip install slideio`.")

        self.path = path
        self.slide = slideio.open_slide(path)

    def read_region(self, location, level, size):
        import numpy as np
        from PIL import Image

        with self.slide.get_scene(0) as scene:
            scaling = 1 / scene.get_zoom_level_info(level).scale
            (x, y) = location
            x_end = min(x + int(scaling * size[0]), scene.size[0])
            y_end = min(y + int(scaling * size[1]), scene.size[1])
            x_width = x_end - x
            y_height = y_end - y
            tile_x = int(np.round(x_width / scaling))
            tile_y = int(np.round(y_height / scaling))
            _tile = scene.read_block((x, y, x_width, y_height), (tile_x, tile_y))
            tile = np.zeros((size[1], size[0], 3), dtype=np.uint8)
            tile[: _tile.shape[0], : _tile.shape[1], :] = _tile[:, :, :3]  # Ensure RGB format
        return Image.fromarray(tile)

    def get_zarr_store(self, tilesize: int = 512):
        from zarr.storage import KVStore

        from ._wsi_store import WsiStore

        return KVStore(WsiStore(self, tilesize))

    def close(self):
        pass

    @property
    def properties(self):
        return {"slideio.objective-power": self.slide.get_scene(0).magnification}

    @property
    def level_downsamples(self):
        return [
            1 / self.slide.get_scene(0).get_zoom_level_info(i).scale
            for i in range(self.slide.get_scene(0).num_zoom_levels)
        ]

    @property
    def level_count(self):
        return self.slide.get_scene(0).num_zoom_levels

    @property
    def level_dimensions(self):
        return [
            (
                self.slide.get_scene(0).get_zoom_level_info(i).size.width,
                self.slide.get_scene(0).get_zoom_level_info(i).size.height,
            )
            for i in range(self.slide.get_scene(0).num_zoom_levels)
        ]

    @property
    def dimensions(self):
        """Get the dimensions of the slide."""
        return self.slide.get_scene(0).size


class XarrayReader(ReaderBase):
    """Reader for xarray data (no backend used)."""

    name = "xarray"

    def __init__(self, slide: DataArray | DataTree):
        self.slide = slide

    def read_region(self, location, level, size):
        """Get a region from the slide."""
        x_start, y_start = location
        x_end, y_end = x_start + size[0], y_start + size[1]
        image = self.slide[f"scale{level}"].image if isinstance(self.slide, DataTree) else self.slide
        tile = image[:, slice(y_start, y_end), slice(x_start, x_end)]
        return tile.transpose("y", "x", "c")

    def close(self):
        """Close the slide."""
        self.slide.close()


def get_reader(backend: str) -> type[ReaderBase]:
    """Get a reader for the specified backend."""
    readers = {
        "openslide": OpenSlideReader,
        "tiffslide": TiffSlideReader,
        "slideio": SlideIOReader,
        "xarray": XarrayReader,
    }
    if backend not in readers:
        raise ValueError(f"Unknown backend: {backend}. Supported backends are: {', '.join(readers.keys())}.")
    return readers[backend]
