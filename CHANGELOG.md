## [1.0.2] - 2024-01-15

### Fix
- When geometries are `GeometryCollection`, convert them back to Polygons (#11)
- Give `min_area` parameter to the right Baysor function in snakemake

### Added
- API tutorial
- `sopa.spatial` tutorial
- Docstrings for the snakemake pipeline utils
- Show right micron scale in the Xenium Explorer

### Changed
- `sopa.stats` is now called `sopa.spatial`

## [1.0.1] - 2024-01-10

### Added

- Tutorial on CLI usage
- Tutorial on image alignment with the Xenium Explorer
- Multi-step segmentation ([#8](https://github.com/gustaveroussy/sopa/issues/8))
- Tutorial for multi-step segmentation and custom segmentation
- Improved installation guide

### Fix
- CLI issue (missing file) when used without Snakemake

### Change

- Use `.parquet` instead of `.zarr.zip` format to store intermediate boundary polygons

## [1.0.0] - 2023-12-26

### Added

- First official release
- Preprint at https://www.biorxiv.org/content/10.1101/2023.12.22.571863v1
