# Spatial-omics pipeline and analysis

<p align="center">
  <img src="./assets/sopa.png" alt="sopa_logo" width="400px"/>
</p>

**S**patial-**o**mics **p**ipeline and **a**nalysis in Python. Built on top of [SpatialData](https://github.com/scverse/spatialdata), it enables processing and analyses of **any** image-based spatial-omics using a standard data structure and output. Sopa was designed for generability and low-memory consumption on large images (scales to `1TB+` images).

The pipeline outputs contain: (i) Xenium Explorer files for interactive visualization, (ii) a HTML report for quick quality controls, and (iii) a SpatialData `.zarr` directory for further analyses.

## Overview

The following illustration describes the main steps of `sopa`:

<p align="center">
  <img src="./assets/overview.png" alt="sopa_overview" width="100%"/>
</p>

## Why use `sopa`

- `sopa` is designed to be memory-efficient, and it scales to large datasets with millions of cells
- Depending on your need, you case use our Snakemake pipeline, our CLI, or our API
- It's very easy to move on to another spatial-omics technology, since `sopa` is general to every image-based spatial-omics
- You can open any data with the [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer), which is a user-friendly software with many functionnalities
- Spatial statistics are optimized, since geometric operations uses `shapely` internally
- You can customize `sopa` and add your own segmentation or annotation tool if desired
