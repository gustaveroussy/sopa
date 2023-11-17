<p align="center">
  <img src="docs/assets/sopa.png" alt="sopa_logo" width="400"/>
</p>

**S**patial-**o**mics **p**ipeline and **a**nalysis in Python. Built on top of [SpatialData](https://github.com/scverse/spatialdata), it enables processing and analyses of **any** imaging-based spatial-omics using a standard data structure and output. Sopa was designed for generability and low-memory consumption on large images (scales to `1TB+` images).

The pipeline outputs are composed of (i) Xenium Explorer files for interactive visualization, (ii) a HTML report for quick quality controls, and (iii) a SpatialData `.zarr` directory for further analyses.

# Documentation

The easiest way to getting started with `sopa` is to check [our documentation](TODO). It contains installation explainations, CLI/API details, and usages examples.

# Overview

The following illustration describes the main steps of `sopa`:

<p align="center">
  <img src="docs/assets/overview.png" alt="sopa_overview" width="100%"/>
</p>

# Installation

### PyPI installation
Sopa can be installed via `PyPI` on all operating system. Make sure you have an environment with `python==3.10`, and run the following command:
```
pip install sopa
```

To install extras (for example if you want to use `cellpose`, or `tangram`), please run:
```
pip install 'sopa[cellpose,baysor,tangram]'
```

Important: even though `pip install 'sopa[baysor]'` will install some dependencies related to baysor, you still have to install the `baysor` command line (see the [official repository](https://github.com/kharchenkolab/Baysor)) if you want to use it.

### Other installation modes

You can clone the repository and run one of these command lines at the root of `sopa`:
```
pip install -e .  # dev mode installation
poetry install    # poetry installation
```

# Features
Sopa comes with three different flavors, each corresponding to a different use case:
- `Snakemake pipeline`: choose a config, and run our pipeline on your spatial data in a couple of minutes
- `CLI`: use our command-line-interface to prototype quickly your own pipeline
- `API`: use directly `sopa` as a python package for full flexibility and customization

### Snakemake pipeline

Clone our repository, choose a config [here](TODO) (or create your own), and execute our pipeline locally or on a high-performance cluster:
```bash
git clone https://github.com/gustaveroussy/sopa.git
cd sopa/workflow
snakemake --configfile=/path/to/yaml_config --config data_path=/path/to/data_directory
```

For more details on `snakemake` configuration and how to properly setup your environments, please refer to the [documentation](TODO).

### CLI

```bash
# Example
> sopa --help # show commands names and arguments
> sopa read merscope_directory --technology merscope # read some data
> sopa patchify merscope_directory.zarr # make patches for low-memory segmentation
> sopa segmentation cellpose merscope_directory.zarr --diameter 60 --channels DAPI # segmentation
> sopa resolve cellpose merscope_directory.zarr # resolve segmentation conflicts at boundaries
> sopa aggregate merscope_directory.zarr --average-intensities # transcripts/channels aggregation
> sopa explorer merscope_directory.zarr # convert for interactive viz
```

For a full description of the API, please refer to the [documentation](TODO).

### API

```python
import sopa

# use the 'sopa' package
```

For a full description of the API, please refer to the [documentation](TODO).

# Cite us
Our article is not published yet. In the meantime, you can cite our **preprint**: TODO