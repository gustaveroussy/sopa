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

## Show the logs

### Show the job/rule logs

Each job has a specific job-id, which is displayed when running `squeue -u $USER`, for instance `2761280`.

By default, the jobs logs are located in `<sopa_repository>/workflow/.snakemake/slurm_logs/rule_<RULE_NAME>/<JOB_ID>.log`.

If the `sopa` repository in located in your workdir (e.g., `/mnt/beegfs/userdata/q_blampey`), you can add this function to your bashrc file:

```sh
function sopa_job_log() {
    find /mnt/beegfs/userdata/$USER/sopa/workflow/.snakemake/slurm_logs -type f -name "$1.log" -exec cat {} \;
}
```

Then, run `sopa_job_log 2761280` to show the logs of the jobs with ID `2761280`.

### Show the main process logs

If you have a snakemake pipeline that is running, but you don't have access to the logs of the main snakemake process, then add the following

```sh
function sopa_smk_log() {
    grep -rl "$1" /mnt/beegfs/userdata/$USER/sopa/workflow/.snakemake/log | xargs cat
}
```

Then, you can copy the job name of any job of the pipeline, for instance `fe76fa78`, and run `sopa_smk_log fe76fa78`.

### Show rule names

You can also look at the rule name of your current running jobs, using `squeue -u $USER -o %i,%P,%.10j,%.40k`.
