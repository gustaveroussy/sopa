To use a Snakemake profile, pass `--workflow-profile profile/my_profile` to your `snakemake` command. Where `my_profile` is an existing profile (currently, `local`, `slurm`, or `lsf`).

In Snakemake>=8.0.0, there are multiple [executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/index.html) that can help you setup a new profile.
