class SopaKeys:
    CELLPOSE_BOUNDARIES = "cellpose_boundaries"
    BAYSOR_BOUNDARIES = "baysor_boundaries"
    PATCHES = "sopa_patches"
    BOUNDS = "bounds"

    UNS_KEY = "sopa_attrs"
    UNS_HAS_TRANSCRIPTS = "transcripts"
    UNS_HAS_INTENSITIES = "intensities"
    UNS_CELL_TYPES = "cell_types"

    INTENSITIES_OBSM = "intensities"

    REGION_KEY = "region"
    SLIDE_KEY = "slide"
    INSTANCE_KEY = "cell_id"
    BAYSOR_DEFAULT_CELL_KEY = "cell"

    Z_SCORES = "z_scores"

    GEOMETRY_AREA = "area"
    GEOMETRY_LENGTH = "length"
    GEOMETRY_ROUNDNESS = "roundness"
    GEOMETRY_COUNT = "n_components"


VALID_DIMENSIONS = ("c", "y", "x")
LOW_AVERAGE_COUNT = 0.01


class ROI:
    KEY = "region_of_interest"
    SCALE_FACTOR = "scale_factor"
    IMAGE_ARRAY_KEY = "image"
    POLYGON_ARRAY_KEY = "polygon"
    IMAGE_KEY = "image_key"
    ELEMENT_TYPE = "element_type"


class SopaFiles:
    SMK_DIR = ".smk_files"
    PATCHES_FILE_CELLPOSE = "patches_file_cellpose"
    PATCHES_DIRS_BAYSOR = "patches_file_baysor"
    BAYSOR_TRANSCRIPTS = "transcripts.csv"
    BAYSOR_CONFIG = "config.toml"
    CELLPOSE_NAME = "cellpose"
    BAYSOR_NAME = "baysor"
