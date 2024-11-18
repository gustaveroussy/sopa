## Installation

### Sopa package

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

    If the Baysor executable is not at `~/.julia/bin/baysor`, please set the `baysor` alias or export the path to the executable via `export baysor=/path/to/exe`.

### Snakemake setup

To use Snakemake, in addition to the above `sopa` environment, you'll need to clone a repository containing the Snakemake workflow:

```sh
git clone https://github.com/gustaveroussy/sopa.git
cd sopa   # move inside the sopa repository
```

Also, make sure you have installed `snakemake`. This does **not** necessarily have to be inside the `sopa` environment: for instance, you can create a new environment specific to snakemake:

```sh
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

Now, follow our [snakemake tutorial](../tutorials/snakemake) to run your first pipeline.

## Usage

Sopa comes in three different flavours, each corresponding to a different use case:

- `Snakemake pipeline`: choose a config, and run our pipeline on your spatial data in a few minutes. See our [snakemake tutorial](../tutorials/snakemake).
- `CLI`: use our [command-line-interface](../tutorials/cli_usage) to prototype quickly your own pipeline
- `API`: use directly `sopa` as a Python package for full flexibility and customization (see a tutorial [here](../tutorials/api_usage))
