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
