On Flamingo (Gustave Roussy Institute), you can use the [slurm-gustave-roussy plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm-gustave-roussy.html) via this config.

## Usage

1. Make sure to have `snakemake>=8.0.0`
2. Install the plugin: `pip install snakemake-executor-plugin-slurm-gustave-roussy`
3. Run Snakemake with the line below (set `data_path` and `configfile` accordingly to your use case):

```sh
snakemake --config data_path=/path/to/data --configfile=config/merscope/base.yaml --workflow-profile profile/slurm-igr
```
