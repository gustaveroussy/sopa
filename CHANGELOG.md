## [1.0.10] - tbd

### Changed
- Default `min_molecules_per_cell` set to 1 in snakemake configs (fix #37)
- Default `min_transcripts` set in snakemake configs

## [1.0.9] - 2024-04-03

### Added:
- Support multiple tables

### Fixed
- Spatial elements not saved when data is not backed

## [1.0.8] - 2024-04-02

Hotfix: resolve issues related to `spatialdata>=1.0.0`

## [1.0.7] - 2024-03-29

### Changed
- Improvements in the CLI and API tutorials
- Sequential segmentation now requires `patchify` to be run independently
- Dependency `spatialdata>=0.1.1`

### Added
- Kwargs can be provided to Cellpose model init

### Fixed
- `set_transformation` issue for image alignment
- Import issue #37 #39

## [1.0.6] - 2024-03-13

### Added
- Spatial join between shapes (`from sdata.spatial import sjoin`)
- H&E tutorial (basic usage)
- New backend for the MERSCOPE reader (requires `rioxarray`, currently experimental, should use less RAM)

### Changed
- Using `MultiscaleSpatialImage` by default for multiplex imaging technologies

### Fixed
- Issue in report creation when channel names are integers

## [1.0.5] - 2024-03-01

### Changed
- Faster image writing for the Xenium Explorer (about x5 speedup)
- Cellpose default model set to `"cyto3"` (new cellpose version)
- Cell GeoDataFrame index consistent with `obs_names`

### Added
- Support for python 3.9 to 3.11 (we still recommend `python==3.10`)
- Support WSI analysis: reader, tissue segmentation, patch embedding (tutorials coming soon)
- Supporting multiple region-of-interest queries
- Can load a custom cellpose model using the `pretrained_model`/`model_type` argument

## [1.0.4] - 2024-02-14

### Fix
- Missing transcript count in cells due to concurrent writing processes (#20)

### Changed
- Explorer images should have a higher contrast (not clipping values anymore)

## [1.0.3] - 2024-02-12

### Breaking changes
- `pixelsize` argument has been renamed to `pixel_size` (the snakemake pipeline is only deprecating it for now)

### Added
- The `phenocycler` reader can now also read `.tif` files (not just `.qptiff`)
- Added missing legend in the HTML report under the "Channels" section (#15)
- The cell area is also stored in the table (in `.obs["area"]`)

### Changed
- The `uniform` toy dataset now has two coordinate systems (better test case)
- Faster table conversion to the Xenium Explorer

### Fixed
- Tight patching more stable with epsilon constant

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
