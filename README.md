<p align="center">
  <img src="docs/assets/sopa.png" alt="sopa_logo" width="400"/>
</p>

**S**patial-**o**mics **p**reprocessing and **a**nalysis in Python. Built on top of [SpatialData](https://github.com/scverse/spatialdata), it provides tools to build a pipeline for any spatial-omics using a standard data structure and output.

# Documentation

The easiest way to getting started with `sopa` is to check [our documentation](TODO). It contains installation explainations, CLI/API details, and usages examples.

# Installation

### PyPI installation
Sopa can be installed via `PyPI` on all operating system. Make sure you have an environment with `python>=3.10`, and run the following command:
```
pip install sopa
```

To install extras (for example if you want to use `cellpose`, or `tangram`), please run:
```
pip install 'sopa[cellpose,tangram]'
```

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

### CLI

### API

# Cite us
TODO