from __future__ import annotations

from pathlib import Path


class WorkflowPaths:
    """
    A class listing the paths to different files that are needed
    or will be created by Snakemake (input and output)
    """

    def __init__(self, config: dict) -> None:
        self.config = config

        ### SpatialData object files
        self.data_path = self.config["data_path"]
        self.sdata_path = Path(self.config["sdata_path"])
        self.sdata_zgroup = self.sdata_path / ".zgroup"  # trick to fix snakemake ChildIOException

        self.shapes_dir = self.sdata_path / "shapes"
        self.points_dir = self.sdata_path / "points"
        self.images_dir = self.sdata_path / "images"
        self.table_dir = self.sdata_path / "tables"

        self.baysor_boundaries = self.shapes_dir / "baysor_boundaries"
        self.cellpose_boundaries = self.shapes_dir / "cellpose_boundaries"

        ### Snakemake cache / intermediate files
        self.sopa_cache = self.sdata_path / ".sopa_cache"

        self.smk_patches = self.sopa_cache / "patches"
        self.smk_aggregation = self.sopa_cache / "aggregation"
        self.smk_table = self.sopa_cache / "table"

        self.smk_cellpose_temp_dir = self.sopa_cache / "cellpose_boundaries"
        self.smk_patches_file_image = self.sopa_cache / "patches_file_image"

        self.smk_transcripts_temp_dir = self.sopa_cache / "transcript_patches"
        self.smk_patches_file_transcripts = self.sopa_cache / "patches_file_transcripts"

        self.smk_cellpose_boundaries = self.sopa_cache / "cellpose_boundaries_done"
        self.smk_baysor_boundaries = self.sopa_cache / "baysor_boundaries_done"
        self.smk_comseg_boundaries = self.sopa_cache / "comseg_boundaries_done"
        self.smk_tissue_segmentation = self.sopa_cache / "tissue_segmentation_done"

        ### Xenium Explorer outputs
        self.explorer_directory = self.sdata_path.with_suffix(".explorer")
        self.explorer_directory.mkdir(parents=True, exist_ok=True)
        self.explorer_experiment = self.explorer_directory / "experiment.xenium"
        self.explorer_image = self.explorer_directory / "morphology.ome.tif"
        self.report = self.explorer_directory / "analysis_summary.html"

        ### Annotation files
        self.annotations = []
        if "annotation" in self.config:
            key = self.config["annotation"].get("args", {}).get("cell_type_key", "cell_type")
            self.annotations = self.table_dir / "table" / "obs" / key

    def temporary_boundaries_paths(self, file_content: str, method_name: str) -> list[str]:
        """Compute the paths to the temporary boundary files

        Args:
            file_content: Content of the file listing the number of patches (or the patches indices)
            method_name: Name of the method: cellpose, baysor, or comseg

        Returns:
            A list of temporary boundary directories or files
        """
        if method_name == "cellpose":
            return [str(self.smk_cellpose_temp_dir / f"{i}.parquet") for i in range(int(file_content))]

        if method_name == "baysor":
            indices = map(int, file_content.split())
            return [str(self.smk_transcripts_temp_dir / str(i) / "segmentation_counts.loom") for i in indices]

        if method_name == "comseg":
            indices = map(int, file_content.split())
            COMSEG_FILES = ["segmentation_polygons.json", "segmentation_counts.h5ad"]
            return [str(self.smk_transcripts_temp_dir / str(i) / file) for i in indices for file in COMSEG_FILES]
