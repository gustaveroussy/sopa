# Flamingo-specific profile

On Flamingo (Gustave Roussy Institute), you can use the [slurm-gustave-roussy plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm-gustave-roussy.html) via this config.

The usage on Flamingo is similar to the [public tutorial](https://gustaveroussy.github.io/sopa/tutorials/snakemake/), but it uses a specific **profile** called `slurm-igr` (see below how to use it).

## Usage

1. Create an environment with `snakemake>=8.0.0`. This can be done with `conda create -c conda-forge -c bioconda -n snakemake snakemake`
2. Activate the environment: `source activate snakemake`
3. Install the plugin: `pip install snakemake-executor-plugin-slurm-gustave-roussy`
4. Move to the `workflow` directory inside the `sopa` repository you cloned.
5. Run Snakemake with the line below (set `data_path` and `configfile` accordingly to your use case):

```sh
snakemake --config data_path=/path/to/data --configfile=config/merscope/base.yaml --workflow-profile profile/slurm-igr
```

## Show job logs

By default, the logs are located in `<sopa_repository>/workflow/.snakemake/slurm_logs/rule_<RULE_NAME>/<JOB_ID>.log`.

If the `sopa` repository in located in your workdir (e.g., `/mnt/beegfs/userdata/q_blampey`), you can add this function to your bashrc file:

```sh
function sopa_log() {
    find /mnt/beegfs/userdata/$USER/sopa/workflow/.snakemake/slurm_logs -type f -name "$1.log" -exec cat {} \;
}
```

Then, run `sopa_log 2761280` to show the logs of the jobs with ID `2761280`.
