rule patch_segmentation_cellpose:
    input:
        paths.smk_patches_file_image,
        paths.smk_patches,
    output:
        paths.temp_dir("cellpose") / "{index}.parquet",
    conda:
        "sopa"
    params:
        cellpose = args["segmentation"]["cellpose"].as_cli(),
        sdata_path = paths.sdata_path,
    shell:
        """
        sopa segmentation cellpose {params.sdata_path} --patch-index {wildcards.index} {params.cellpose}
        """


rule resolve_cellpose:
    input:
        get_input_resolve("image", "cellpose"),
    output:
        touch(paths.segmentation_done("cellpose")),
    conda:
        "sopa"
    params:
        sdata_path = paths.sdata_path,
    shell:
        """
        sopa resolve cellpose {params.sdata_path}
        """
