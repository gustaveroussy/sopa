def get_reader(backend: str):
    """Get a reader for the specified backend."""
    if backend == "openslide":
        return OpenSlideReader
    elif backend == "tiffslide":
        return TiffSlideReader
    elif backend == "slideio":
        return SlideIOReader
    else:
        raise ValueError(
            f"Unsupported backend: {backend}. Supported backends are 'openslide', 'tiffslide', and 'slideio'."
        )


class ReaderBase:
    path: str
    name: "base"
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

    def __init__(self, path: str):
        from openslide import OpenSlide

        self.name = "openslide"
        self.path = path
        self.slide = OpenSlide(path)

    def read_region(self, location, level, size):
        """Get a region from the slide."""
        return self.slide.read_region(location, level, size)

    def get_zarr_store(self, tilesize: int = 512):
        from ._wsi_store import WsiStore
        from zarr.storage import KVStore

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

    def __init__(self, path: str):
        import tiffslide

        self.name = "tiffslide"
        self.path = path
        self.slide = tiffslide.open_slide(path)

    def read_region(self, location, level, size):
        """Get a region from the slide."""
        return self.slide.read_region(location, level, size)

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


    def __init__(self, path: str):
        import slideio

        self.name = "slideio"
        self.path = path
        self.slide = slideio.open_slide(path)

    def read_region(self, location, level, size):
        import numpy as np
        
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
            tile = np.zeros((size[0], size[1], scene.num_channels), dtype=np.uint8)
            tile[: _tile.shape[0], : _tile.shape[1], :] = _tile
        return np.array(tile)

    def get_zarr_store(self, tilesize: int = 512):
        from ._wsi_store import WsiStore
        from zarr.storage import KVStore

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
