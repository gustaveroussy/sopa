class SopaKeys:
    # Segmentation keys
    CELLPOSE_BOUNDARIES = "cellpose_boundaries"
    BAYSOR_BOUNDARIES = "baysor_boundaries"
    COMSEG_BOUNDARIES = "comseg_boundaries"
    CUSTOM_BOUNDARIES = "custom_boundaries"

    # Patches keys
    PATCHES = "image_patches"
    TRANSCRIPTS_PATCHES = "transcripts_patches"
    EMBEDDINGS_PATCHES = "embeddings_patches"
    CACHE_PATH_KEY = "cache_path"
    BOUNDS = "bboxes"
    PATCHES_ILOCS = "ilocs"
    ROI = "region_of_interest"
    PRIOR_SHAPES_KEY = "prior_shapes_key"
    POINTS_KEY = "points_key"

    # Other SpatialData keys
    TABLE = "table"
    OLD_TABLE_PREFFIX = "old_"

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
    BOUNDARIES = "boundaries_shapes"
    GENE_EXCLUDE_PATTERN = "gene_exclude_pattern"  # usage: columns.str.match(pattern, case=False)
    UID = "sopa_uid"


class SopaFiles:
    SOPA_CACHE_DIR = ".sopa_cache"
    TRANSCRIPT_CACHE_DIR = "transcript_patches"
    PATCHES_FILE_IMAGE = "patches_file_image"
    PATCHES_FILE_TRANSCRIPTS = "patches_file_transcripts"
    TRANSCRIPTS_FILE = "transcripts.csv"
    CENTROIDS_FILE = "centroids.csv"
    JSON_CONFIG_FILE = "config.json"
    TOML_CONFIG_FILE = "config.toml"


VALID_DIMENSIONS = ("c", "y", "x")
LOW_AVERAGE_COUNT = 0.01
EPS = 1e-5
ATTRS_KEY = "spatialdata_attrs"
