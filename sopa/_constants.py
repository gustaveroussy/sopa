class SopaKeys:
    CELLPOSE_BOUNDARIES = "cellpose_boundaries"
    BAYSOR_BOUNDARIES = "baysor_boundaries"
    PATCHES = "sopa_patches"
    BOUNDS = "bounds"

    INTENSITIES_OBSM = "intensities"

    REGION_KEY = "region"
    SLIDE_KEY = "slide"
    INSTANCE_KEY = "cell_id"
    BAYSOR_DEFAULT_CELL_KEY = "cell"

    Z_SCORES = "z_scores"


VALID_DIMENSIONS = ("c", "y", "x")


class ROI:
    KEY = "region_of_interest"
    SCALE_FACTOR = "scale_factor"
    IMAGE_ARRAY_KEY = "image"
    POLYGON_ARRAY_KEY = "polygon"
    IMAGE_KEY = "image_key"
    ELEMENT_TYPE = "element_type"


class SopaFiles:
    SMK_DIR = ".smk_files"
    NUM_PATCHES_CELLPOSE = "n_patches_cellpose"
    PATCHES_DIRS_BAYSOR = "n_patches_baysor"
    BAYSOR_TRANSCRIPTS = "transcripts.csv"
    BAYSOR_CONFIG = "config.toml"
    CELLPOSE_NAME = "cellpose"
    BAYSOR_NAME = "baysor"
