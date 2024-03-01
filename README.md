<p align="center">
  <img src="https://raw.githubusercontent.com/gustaveroussy/sopa/master/docs/assets/sopa.png" alt="sopa_logo" width="250"/>
</p>

# Spatial-omics pipeline and analysis
[![PyPI](https://img.shields.io/pypi/v/sopa.svg)](https://pypi.org/project/sopa)
[![Downloads](https://static.pepy.tech/badge/sopa)](https://pepy.tech/project/sopa)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://gustaveroussy.github.io/sopa)
![Build](https://github.com/gustaveroussy/sopa/workflows/ci/badge.svg)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)
[![License](https://img.shields.io/pypi/l/sopa.svg)](https://github.com/gustaveroussy/sopa/blob/master/LICENSE)
[![Imports: isort](https://img.shields.io/badge/imports-isort-blueviolet)](https://pycqa.github.io/isort/)

Built on top of [SpatialData](https://github.com/scverse/spatialdata), Sopa enables processing and analyses of image-based spatial-omics using a standard data structure and output. We currently support the following technologies: Xenium, MERSCOPE, CosMX, PhenoCycler, MACSima, Hyperion. Sopa was designed for generability and low memory consumption on large images (scales to `1TB+` images).

The pipeline outputs contain: (i) Xenium Explorer files for interactive visualization, (ii) an HTML report for quick quality controls, and (iii) a SpatialData `.zarr` directory for further analyses.

# Documentation

The easiest way to start with `sopa` is to check [our documentation](https://gustaveroussy.github.io/sopa). It contains installation explanations, CLI/API details, and tutorials.

# Overview

The following illustration describes the main steps of `sopa`:

<p align="center">
  <img src="https://raw.githubusercontent.com/gustaveroussy/sopa/master/docs/assets/overview_white.png" alt="sopa_overview" width="100%"/>
</p>

# Installation

### PyPI installation
Sopa can be installed via `PyPI` on all operating systems. The preferred Python version is `python==3.10`, but we also support `3.9` to `3.11`. On a new environment, run the following command:
```
pip install sopa
```

To install extras (for example, if you want to use `snakemake`/`cellpose`/`baysor`/`tangram`), please run:
```
pip install 'sopa[snakemake,cellpose,baysor,tangram]'
```

Important: even though `pip install 'sopa[baysor]'` will install some dependencies related to baysor, you still have to install the `baysor` command line (see the [official repository](https://github.com/kharchenkolab/Baysor)) if you want to use it.

### Other installation modes

You can clone the repository and run one of these command lines at the root of `sopa`:
```
pip install -e . # dev mode installation
poetry install    # poetry installation
```

# Features
Sopa comes in three different flavours, each corresponding to a different use case:
- `Snakemake pipeline`: choose a config, and run our pipeline on your spatial data in a couple of minutes
- `CLI`: use our command-line-interface for prototyping quickly your own pipeline
- `API`: use directly `sopa` as a Python package for complete flexibility and customization

### Snakemake pipeline

Clone our repository, choose a config [here](https://github.com/gustaveroussy/sopa/tree/master/workflow/config) (or create your own), and execute our pipeline locally or on a high-performance cluster:
```bash
git clone https://github.com/gustaveroussy/sopa.git
cd sopa/workflow
snakemake --configfile=/path/to/yaml_config --config data_path=/path/to/data_directory --cores 1 --use-conda
```

For more details on `snakemake` configuration and how to properly setup your environments, please refer to the [documentation](https://gustaveroussy.github.io/sopa/tutorials/snakemake/).

### CLI

Below are examples of commands that can be run with the `sopa` CLI:

```bash
> sopa --help # show command names and arguments
> sopa read merscope_directory --technology merscope # read some data
> sopa patchify image merscope_directory.zarr # make patches for low-memory segmentation
> sopa segmentation cellpose merscope_directory.zarr --diameter 60 --channels DAPI # segmentation
> sopa resolve cellpose merscope_directory.zarr # resolve segmentation conflicts at boundaries
> sopa aggregate merscope_directory.zarr --average-intensities # transcripts/channels aggregation
> sopa explorer write merscope_directory.zarr # convert for interactive vizualisation
```

For a complete description of the CLI, please refer to the [documentation](https://gustaveroussy.github.io/sopa/cli).

### API

```python
import sopa

# use the 'sopa' python package
```

For a complete API description, please refer to the [documentation](https://gustaveroussy.github.io/sopa).

# Cite us
Our article is not published yet. In the meantime, you can cite our [preprint](https://www.biorxiv.org/content/10.1101/2023.12.22.571863v1):

```txt
@article {Blampey2023.12.22.571863,
    author = {Quentin Blampey & Kevin Mulder et al.},
    title = {Sopa: a technology-invariant pipeline for analyses of image-based spatial-omics},
    elocation-id = {2023.12.22.571863},
    year = {2023},
    doi = {10.1101/2023.12.22.571863},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2023/12/23/2023.12.22.571863},
    eprint = {https://www.biorxiv.org/content/early/2023/12/23/2023.12.22.571863.full.pdf},
    journal = {bioRxiv}
}
```
