<p align="center">
  <img src="https://raw.githubusercontent.com/gustaveroussy/sopa/main/docs/assets/sopa.png" alt="sopa_logo" width="250"/>
</p>
<p align="center"><b><i>
	Spatial omics pipeline and analysis
</b></i></p>

<div align="center">

[![PyPI](https://img.shields.io/pypi/v/sopa.svg)](https://pypi.org/project/sopa)
[![Downloads](https://static.pepy.tech/badge/sopa)](https://pepy.tech/project/sopa)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://gustaveroussy.github.io/sopa)
![Build](https://github.com/gustaveroussy/sopa/workflows/ci/badge.svg)
[![License](https://img.shields.io/pypi/l/sopa.svg)](https://github.com/gustaveroussy/sopa/blob/main/LICENSE)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

</div>

Built on top of [SpatialData](https://github.com/scverse/spatialdata), Sopa enables processing and analyses of spatial omics data with single-cell resolution (spatial transcriptomics or multiplex imaging data) using a standard data structure and output. We currently support the following technologies: Xenium, Visium HD, MERSCOPE, CosMX, PhenoCycler, MACSima, Molecural Cartography, and others. Sopa was designed for generability and low memory consumption on large images (scales to `1TB+` images).

🎉 `sopa==2.0.0` is out! Check [our migration guide](https://github.com/gustaveroussy/sopa/discussions/138) to smoothly update your code base.

## Documentation

Check [Sopa's documentation](https://gustaveroussy.github.io/sopa) to get started. It contains installation explanations, CLI/API details, and tutorials.

## Overview

The following illustration describes the main steps of `sopa`:

<p align="center">
  <img src="https://raw.githubusercontent.com/gustaveroussy/sopa/main/docs/assets/overview_white.png" alt="sopa_overview" width="100%"/>
</p>

## Installation

Sopa can be installed via `PyPI` on all operating systems, with the only requirement being Python (`>=3.10` and `<3.13`). On a new environment, run the following command:
```sh
pip install sopa
```

See the [installation section](https://gustaveroussy.github.io/sopa/getting_started/) from the docs for more details, to install extras or to use other installations modes.

## Features
Sopa comes in three different flavours, each corresponding to a different use case:
- `API`: use directly `sopa` as a Python package for complete flexibility and customization
- `Snakemake pipeline`: choose a config, and run our pipeline on your spatial data in a couple of minutes
- `CLI`: use our command-line-interface for prototyping quickly your own pipeline

### API

Below is an example of a minimal API usage. For a complete API description, please refer to the [documentation](https://gustaveroussy.github.io/sopa).

```python
import sopa

sdata = sopa.io.xenium("path/to/data") # reading Xenium data

sopa.make_image_patches(sdata) # creating overlapping patches
sopa.segmentation.cellpose(sdata, "DAPI", diameter=30) # running cellpose segmentation
sopa.aggregate(sdata) # counting the transcripts inside the cells
```

### Snakemake pipeline

Clone our repository, choose a config [here](https://github.com/gustaveroussy/sopa/tree/main/workflow/config) (or create your own), and execute our pipeline locally or on a high-performance cluster:
```bash
git clone https://github.com/gustaveroussy/sopa.git
cd sopa/workflow
snakemake --configfile=/path/to/yaml_config --config data_path=/path/to/data_directory --cores 1 --use-conda
```

For more details on `snakemake` configuration and how to properly setup your environments, please refer to the [documentation](https://gustaveroussy.github.io/sopa/tutorials/snakemake/).

### CLI

Below are examples of commands that can be run with the `sopa` CLI. For a complete description of the CLI, please refer to the [documentation](https://gustaveroussy.github.io/sopa/cli).

```bash
> sopa --help # show command names and arguments
> sopa convert merscope_directory --technology merscope # read some data
> sopa patchify image merscope_directory.zarr # make patches for low-memory segmentation
> sopa segmentation cellpose merscope_directory.zarr --diameter 60 --channels DAPI # segmentation
> sopa resolve cellpose merscope_directory.zarr # resolve segmentation conflicts at boundaries
> sopa aggregate merscope_directory.zarr --average-intensities # transcripts/channels aggregation
> sopa explorer write merscope_directory.zarr # convert for interactive vizualisation
```

## Cite us
Our article is published in [Nature Communications](https://www.nature.com/articles/s41467-024-48981-z). You can cite Sopa as below:

```txt
@article{blampey_sopa_2024,
	title = {Sopa: a technology-invariant pipeline for analyses of image-based spatial omics},
	volume = {15},
	url = {https://www.nature.com/articles/s41467-024-48981-z},
	doi = {10.1038/s41467-024-48981-z},
	journal = {Nature Communications},
	author = {Blampey, Quentin and Mulder, Kevin and Gardet, Margaux and Christodoulidis, Stergios and Dutertre, Charles-Antoine and André, Fabrice and Ginhoux, Florent and Cournède, Paul-Henry},
	year = {2024},
	note = {Publisher: Nature Publishing Group},
	pages = {4981},
}
```
