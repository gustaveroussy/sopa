class SopaKeys:
    # Segmentation keys
    CELLPOSE_BOUNDARIES = "cellpose_boundaries"
    BAYSOR_BOUNDARIES = "baysor_boundaries"
    COMSEG_BOUNDARIES = "comseg_boundaries"

    # Patches keys
    PATCHES = "sopa_patches"
    TRANSCRIPT_PATCHES = "sopa_patches_transcripts"
    PATCHES_INFERENCE_KEY = "sopa_patches_inference"
    CACHE_PATH_KEY = "cache_path"
    BOUNDS = "bboxes"
    PATCHES_ILOCS = "ilocs"

    # Other SpatialData keys
    TABLE = "table"
    OLD_TABLE = "old_table"

    # Table keys
    UNS_KEY = "sopa_attrs"
    UNS_HAS_TRANSCRIPTS = "transcripts"
    UNS_HAS_INTENSITIES = "intensities"
    UNS_CELL_TYPES = "cell_types"
    INTENSITIES_OBSM = "intensities"
    REGION_KEY = "region"
    SLIDE_KEY = "slide"
    INSTANCE_KEY = "cell_id"
    DEFAULT_CELL_KEY = "cell"
    ORIGINAL_AREA_OBS = "baysor_area"
    CELL_OVERLAY_KEY = "is_overlay"
    AREA_OBS = "area"
    Z_SCORES = "z_scores"

    # Geometry keys
    GEOMETRY_AREA = "area"
    GEOMETRY_LENGTH = "length"
    GEOMETRY_ROUNDNESS = "roundness"
    GEOMETRY_COUNT = "n_components"


class SopaAttrs:
    CELL_SEGMENTATION = "cell_segmentation_image"
    TISSUE_SEGMENTATION = "tissue_segmentation_image"
    BINS_TABLE = "bins_table"
    TRANSCRIPTS = "transcripts_dataframe"
    GENE_COLUMN = "feature_key"


class ROI:
    KEY = "region_of_interest"
    SCALE_FACTOR = "scale_factor"
    IMAGE_ARRAY_KEY = "image"
    POLYGON_ARRAY_KEY = "polygon"
    IMAGE_KEY = "image_key"
    ELEMENT_TYPE = "element_type"


class SopaFiles:
    SOPA_CACHE_DIR = ".sopa_cache"
    TRANSCRIPT_TEMP_DIR = "transcript_patches"
    PATCHES_FILE_IMAGE = "patches_file_image"
    PATCHES_DIRS_BAYSOR = "patches_file_baysor"
    PATCHES_DIRS_COMSEG = "patches_file_comseg"
    TRANSCRIPTS_FILE = "transcripts.csv"
    CENTROIDS_FILE = "centroids.csv"
    JSON_CONFIG_FILE = "config.json"
    TOML_CONFIG_FILE = "config.toml"


VALID_DIMENSIONS = ("c", "y", "x")
LOW_AVERAGE_COUNT = 0.01
EPS = 1e-5
