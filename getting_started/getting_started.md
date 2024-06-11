## Installation

Sopa can be installed on every OS with `pip` or [`poetry`](https://python-poetry.org/docs/).

The preferred Python version is `python==3.10`, but we also support `3.9` to `3.11`.

!!! note "Advice (optional)"

    We advise creating a new environment via a package manager (except if you use Poetry, which will automatically create the environment).

    For instance, you can create a new `conda` environment:

    ```bash
    conda create --name sopa python=3.10
    conda activate sopa
    ```

Choose one of the following, depending on your needs (it should take at most a few minutes):

=== "From PyPI"

    ``` bash
    pip install sopa

    # or to install extras
    pip install 'sopa[cellpose,baysor,tangram]'
    ```

=== "Local install (pip)"

    ``` bash
    git clone https://github.com/gustaveroussy/sopa.git
    cd sopa

    pip install .
    ```

=== "Poetry (dev mode)"

    ``` bash
    git clone https://github.com/gustaveroussy/sopa.git
    cd sopa

    poetry install --all-extras
    ```

!!! warning "Baysor usage"
    Even though `pip install 'sopa[baysor]'` will install some dependencies related to baysor, you still have to install the `baysor` command line (see the [official repository](https://github.com/kharchenkolab/Baysor)) if you want to use it.

## Snakemake setup

To use the Snakemake pipeline, the installation process is slightly different because you'll need the whole repository.

1. Clone the `sopa` repository, and move to the root of the project:
```sh
git clone https://github.com/gustaveroussy/sopa.git
cd sopa
```

1. Create a `conda` environment called `sopa`:
```sh
conda create --name sopa python=3.10
```

1. At the root of the `sopa` directory, install the package in dev mode, and choose the extras you want (among cellpose/baysor/tangram, depending on your desired usage):
```sh
conda activate sopa
pip install -e ".[snakemake,cellpose,baysor,tangram]"
```

Now, follow our [snakemake tutorial](../tutorials/snakemake) to run your first pipeline.

!!! Note
    You can also use a separate environment for `snakemake`. In this case, you don't need to install the `'snakemake'` extra when installing `sopa`. But you may still need to install other extras, for instance, `'cellpose'` if you plan to run Cellpose.

## Usage

Sopa comes in three different flavours, each corresponding to a different use case:

- `Snakemake pipeline`: choose a config, and run our pipeline on your spatial data in a few minutes. See our [snakemake tutorial](../tutorials/snakemake).
- `CLI`: use our [command-line-interface](../tutorials/cli_usage) to prototype quickly your own pipeline
- `API`: use directly `sopa` as a Python package for full flexibility and customization (see a tutorial [here](../tutorials/api_usage))
