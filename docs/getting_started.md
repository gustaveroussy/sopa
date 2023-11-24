## Installation

Sopa can be installed on every OS with `pip` or [`poetry`](https://python-poetry.org/docs/).

For now, `python==3.10` is required, but more versions will be supported soon.

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
    pip install 'sopa[snakemake,cellpose,baysor,tangram]'
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

## Usage

Sopa comes with three different flavors, each corresponding to a different use case:

- `Snakemake pipeline`: choose a config, and run our pipeline on your spatial data in a couple of minutes. See our [snakemake guide](../pipeline).
- `CLI`: use our [command-line-interface](../cli) to prototype quickly your own pipeline
- `API`: use directly `sopa` as a python package for full flexibility and customization (see the API documentation on the sidebar)