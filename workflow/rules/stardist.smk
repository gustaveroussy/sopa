rule patch_segmentation_stardist:
    input:
        paths.smk_patches_file_image,
        paths.smk_patches,
    output:
        paths.temp_dir("stardist") / "{index}.parquet",
    conda:
        "sopa"
    params:
        stardist = args["segmentation"]["stardist"].as_cli(),
        sdata_path = paths.sdata_path,
    shell:
        """
        sopa segmentation stardist {params.sdata_path} --patch-index {wildcards.index} {params.stardist}
        """


rule resolve_stardist:
    input:
        get_input_resolve("image", "stardist"),
    output:
        touch(paths.segmentation_done("stardist")),
    conda:
        "sopa"
    params:
        sdata_path = paths.sdata_path,
    shell:
        """
        sopa resolve stardist {params.sdata_path}
        """
