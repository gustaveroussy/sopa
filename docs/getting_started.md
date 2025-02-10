## Installation

### Sopa package

Sopa can be installed on every OS with `pip` or [`poetry`](https://python-poetry.org/docs/).

The preferred Python version is `python==3.10`, but we also support `3.11` and `3.12`.

!!! note "Advice (optional)"

    We advise creating a new environment via a package manager (except if you use Poetry, which will automatically create the environment).

    For instance, you can create a new `conda` environment:

    ```bash
    conda create --name sopa python=3.10
    conda activate sopa
    ```

Choose one of the following, depending on your needs (it should take at most a few minutes):

=== "From PyPI"

    ```sh
    pip install sopa
    ```

    To install extras (for example, if you want to use `cellpose`/`baysor`), please run:

    ```sh
    # choose any valid extra among cellpose/baysor/tangram/wsi
    pip install 'sopa[cellpose,baysor]'
    ```

=== "Editable mode"

    ``` bash
    git clone https://github.com/gustaveroussy/sopa.git
    cd sopa

    # no extra
    pip install  -e .

    # or, to install extras, among cellpose/baysor/tangram/wsi:
    pip install -e '.[cellpose,baysor]'
    ```

=== "Poetry (dev mode)"

    ``` bash
    git clone https://github.com/gustaveroussy/sopa.git
    cd sopa

    poetry install --all-extras
    ```

!!! warning "Baysor usage"
    Even though `pip install 'sopa[baysor]'` will install some dependencies related to baysor, you still have to install the `baysor` command line (see the [official documentation](https://kharchenkolab.github.io/Baysor/dev/installation/)) if you want to use it.

    If the Baysor executable is not at `~/.julia/bin/baysor`, please make the `baysor` command available (e.g., via creating a symlink `~/.local/bin/baysor` pointing to executable), or export the path to the executable via `export baysor=/path/to/baysor/executable`.

### Snakemake setup

To use Snakemake, in addition to the above `sopa` environment, you'll need to clone a repository containing the Snakemake workflow:

```sh
git clone https://github.com/gustaveroussy/sopa.git
cd sopa   # move inside the sopa repository
```

Also, make sure you have installed `snakemake>=8.0.0`. This does **not** necessarily have to be inside the `sopa` environment: for instance, you can create a new environment specific to snakemake:

```sh
# this will create a new environment called "snakemake"
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

Now, follow our [snakemake tutorial](../tutorials/snakemake) to run your first pipeline.

## Usage

Sopa comes in three different flavours, each corresponding to a different use case:

- `API`: use directly `sopa` as a Python package for full flexibility and customization (see a tutorial [here](../tutorials/api_usage))
- `Snakemake pipeline`: choose a config, and run our pipeline on your spatial data in a few minutes. See our [snakemake tutorial](../tutorials/snakemake).
- `CLI`: use our [command-line-interface](../tutorials/cli_usage) to prototype quickly your own pipeline
