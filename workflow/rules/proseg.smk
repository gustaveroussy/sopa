rule patch_segmentation_proseg:
    input:
        paths.smk_patches_file_transcripts,
    output:
        touch(paths.segmentation_done("proseg")),
        touch(paths.smk_table),
    conda:
        "sopa"
    params:
        proseg_config = args["segmentation"]["proseg"].as_cli(keys=["command_line_suffix"]),
        sdata_path = paths.sdata_path,
    shell:
        """
        sopa segmentation proseg {params.sdata_path} {params.proseg_config}
        """
