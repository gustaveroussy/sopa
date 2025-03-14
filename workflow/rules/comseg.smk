rule patch_segmentation_comseg:
    input:
        paths.smk_patches_file_transcripts,
    output:
        paths.smk_transcripts_temp_dir / "{index}" / "segmentation_polygons.json",
        paths.smk_transcripts_temp_dir / "{index}" / "segmentation_counts.h5ad",
    conda:
        "sopa"
    params:
        comseg_config = args["segmentation"]["comseg"].as_cli(keys=["config"]),
        sdata_path = paths.sdata_path,
    shell:
        """
        sopa segmentation comseg {params.sdata_path} --patch-index {wildcards.index} {params.comseg_config}
        """


rule resolve_comseg:
    input:
        files = get_input_resolve("transcripts", "comseg"),
    output:
        touch(paths.segmentation_done("comseg")),
        touch(paths.smk_table),
    conda:
        "sopa"
    params:
        resolve = args.resolve_transcripts(),
        sdata_path = paths.sdata_path,
        smk_transcripts_temp_dir = paths.smk_transcripts_temp_dir,
    shell:
        """
        sopa resolve comseg {params.sdata_path} {params.resolve}

        rm -r {params.smk_transcripts_temp_dir}    # cleanup large comseg files
        """
