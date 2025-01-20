# Spatial-omics pipeline and analysis

<p align="center">
  <img src="./assets/sopa.png" alt="sopa_logo" width="250px"/>
</p>

Built on top of [SpatialData](https://github.com/scverse/spatialdata), Sopa enables processing and analyses of spatial omics data with single-cell resolution (spatial transcriptomics or multiplex imaging data) using a standard data structure and output. We currently support the following technologies: Xenium, Visium HD, MERSCOPE, CosMX, PhenoCycler, MACSima, Hyperion. Sopa was designed for generability and low memory consumption on large images (scales to `1TB+` images).

🎉 `sopa==2.0.0` is out! It introduces many new API features; check [our migration guide](https://github.com/gustaveroussy/sopa/discussions/138) to smoothly update your code base.

## Overview

The following illustration describes the main steps of `sopa`:

<p align="center">
  <img src="./assets/overview.png" alt="sopa_overview" width="100%"/>
</p>

## Why use `sopa`

- `sopa` is designed to be memory-efficient, and it scales to large datasets with millions of cells
- It's straightforward to move on to another spatial omics technology since `sopa` is general to every spatial omics with single-cell resolution
- Depending on your need, you can use our Snakemake pipeline, our CLI, or our API
- You can open any data with the [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer), which is a user-friendly software with many functions
- Spatial operations are optimized since geometric operations use `shapely` internally
- You can customize `sopa` and add your own segmentation or annotation tool if desired
- `sopa` integrates naturally with other community tools such as [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) or [Squidpy](https://squidpy.readthedocs.io/en/latest/index.html).
