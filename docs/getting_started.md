## Installing Sopa

Sopa can be installed on every OS with `pip`, from `python>=3.10` to `python<3.13`.

!!! note "Advice (optional)"

    We advise creating a new environment via a package manager (except if you use Poetry, which will automatically create the environment).

    For instance, you can create a new `conda` environment:

    ```bash
    conda create --name sopa python=3.10
    conda activate sopa
    ```

Choose one of the following, depending on your needs:

=== "From PyPI"

    ```sh
    pip install sopa
    ```

=== "Editable mode"

    ``` bash
    git clone https://github.com/gustaveroussy/sopa.git
    cd sopa

    # no extra
    pip install  -e .

    # or, to install extras, among cellpose/baysor/stardist/wsi:
    pip install -e '.[cellpose,baysor]'
    ```

=== "uv (dev mode)"

    ``` bash
    git clone https://github.com/gustaveroussy/sopa.git
    cd sopa

    uv sync --all-extras --dev
    ```

!!! warning "Extra dependencies"
    Dependending on the segmentation tool that you use, you'll need extras, as detailed in the next section.

## Extra dependencies

By default, `sopa` only install the minimal dependencies to avoid a heavy installation. Depending on your usage, you can install some extras. The available extras are listed below.

=== "Cellpose"
    If you need to run Cellpose, you can use the corresponding extra:

    ```sh
    pip install 'sopa[cellpose]'

    # you can also combine extras: pip install 'sopa[cellpose,baysor,wsi,stardist]'
    ```

=== "Baysor"
    To use [Baysor](https://kharchenkolab.github.io/Baysor/dev/), you'll first need to install Sopa with the `baysor` extra:

    ```sh
    pip install 'sopa[baysor]'
    ```

    **Important**: then, you also have to install the `baysor` command line as detailed in the [official documentation](https://kharchenkolab.github.io/Baysor/dev/installation/).

    !!! info "Executable path"
        If the Baysor executable is not at `~/.julia/bin/baysor`, please make the `baysor` command available (e.g., via creating a symlink `~/.local/bin/baysor` pointing to the executable), or export the path to the executable via `export baysor=/path/to/baysor/executable`.


=== "Proseg"
    [Proseg](https://github.com/dcjones/proseg) has to be installed independently, this can be done with [`cargo`](https://doc.rust-lang.org/cargo/getting-started/installation.html):

    ```sh
    cargo install proseg
    ```

    !!! info "Executable path"
        If the `proseg` executable is not at `~/.cargo/bin/proseg`, please make the `proseg` command available (e.g., via creating a symlink `~/.local/bin/proseg` pointing to the executable), or export the path to the executable via `export proseg=/path/to/proseg/executable`.
=== "Stardist"
    If you need to run [Stardist](https://github.com/stardist/stardist), you can install the corresponding extra:

    ```sh
    pip install 'sopa[stardist]'
    ```
=== "Comseg"
    If you need to run Comseg, you can install it via pip:

    ```sh
    pip install comseg
    ```
=== "WSI"
    If you need to work on whole slide images / H&E images, you can install the corresponding extras as below. You can also consider the `stardist` extra, if you want to run cell segmentation on the H&E image.

    ```sh
    pip install 'sopa[wsi]'
    ```

## Snakemake setup

If you plan to use Snakemake, in addition to the above `sopa` environment, you'll need to clone the `sopa` repository containing the Snakemake workflow:

```sh
git clone https://github.com/gustaveroussy/sopa.git
cd sopa   # move inside the sopa repository
```

Also, make sure you have installed `snakemake>=8.0.0`. This does **not** necessarily have to be inside the `sopa` environment — for instance, you can create a new environment specific to snakemake:

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
