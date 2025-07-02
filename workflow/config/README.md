# Pipeline config files

Sopa's config are `.yaml` files, they work for both `snakemake` and `nextflow`. We provide many examples of such a config in this directory.
You can choose an existing one (inside the directory corresponding to your technology of interest) or make your own config.

All arguments are detailed in the [`example_commented.yaml`](https://github.com/gustaveroussy/sopa/blob/main/workflow/config/example_commented.yaml) file.

### Call the online config files from Nextflow

If running Sopa on `nextflow`, you can directly use the config files from GitHub without downloading them.

For instance, if you want to use the Xenium proseg config at `https://github.com/gustaveroussy/sopa/blob/main/workflow/config/xenium/proseg.yaml`, just click on the "Raw" button from the GitHub interface (top right) to access its raw content.

Then, pass this "raw-content" URL to the `--configfile` argument when running `nextflow`, for instance:

```sh
--configfile https://raw.githubusercontent.com/gustaveroussy/sopa/refs/heads/main/workflow/config/xenium/proseg.yaml
```
