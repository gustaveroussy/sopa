rule patch_segmentation_proseg:
    input:
        (self.segmentation_done("stardist") if args.use("stardist") else []) if config["read"]["technology"] == "visium_hd" else paths.smk_patches_file_transcripts,
        paths.sdata_zgroup,
        paths.segmentation_done("tissue") if args.use("tissue") else [],
    output:
        touch(paths.segmentation_done("proseg")),
        touch(paths.smk_table),
    conda:
        "sopa"
    params:
        proseg_config = args["segmentation"]["proseg"].as_cli(exclude=[] if config["read"]["technology"] == "visium_hd" else ["prior_shapes_key"]),
        sdata_path = paths.sdata_path,
    shell:
        """
        sopa segmentation proseg {params.sdata_path} {params.proseg_config}
        """
