rule patch_segmentation_baysor:
    input:
        paths.smk_patches_file_transcripts,
    output:
        paths.smk_transcripts_temp_dir / "{index}" / "segmentation_counts.loom",
    conda:
        "sopa"
    params:
        baysor_config = args["segmentation"]["baysor"].as_cli(keys=["config"]),
        sdata_path = paths.sdata_path,
    resources:
        cpus_per_task=1,
    shell:
        """
        export JULIA_NUM_THREADS={resources.cpus_per_task} # parallelize within each patch for Baysor >= v0.7


        if command -v module &> /dev/null; then
            module purge

            if module avail baysor &> /dev/null; then
                module load baysor
            fi
        fi

        sopa segmentation baysor {params.sdata_path} --patch-index {wildcards.index} {params.baysor_config}
        """

rule resolve_baysor:
    input:
        files = get_input_resolve("transcripts", "baysor"),
    output:
        touch(paths.segmentation_done("baysor")),
        touch(paths.smk_table),
    conda:
        "sopa"
    params:
        resolve = args.resolve_transcripts(),
        sdata_path = paths.sdata_path,
        smk_transcripts_temp_dir = paths.smk_transcripts_temp_dir,
    shell:
        """
        sopa resolve baysor {params.sdata_path} {params.resolve}

        rm -r {params.smk_transcripts_temp_dir}    # cleanup large baysor files
        """
