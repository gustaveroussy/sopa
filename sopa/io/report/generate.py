import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import spatialdata
from spatialdata import SpatialData

from ..._constants import LOW_AVERAGE_COUNT, SopaKeys
from ..._sdata import (
    get_boundaries,
    get_intensities,
    get_intrinsic_cs,
    get_spatial_image,
)
from .engine import (
    CodeBlock,
    Columns,
    Image,
    Message,
    Paragraph,
    ProgressBar,
    Root,
    Section,
    SubSection,
)


def write_report(path: str, sdata: SpatialData):
    sections = SectionBuilder(sdata).compute_sections()

    Root(sections).write(path)


class SectionBuilder:
    def __init__(self, sdata: SpatialData):
        self.sdata = sdata

    def _table_has(self, key, default=False):
        if not SopaKeys.UNS_KEY in self.sdata.table.uns:
            return default
        return self.sdata.table.uns[SopaKeys.UNS_KEY].get(key, default)

    def general_section(self):
        return Section(
            "General",
            [
                SubSection(
                    "SpatialData information",
                    [
                        Paragraph(
                            f"Sopa is using <a href='https://spatialdata.scverse.org/en/latest/'>SpatialData</a> under the hood. This is how the object looks like:"
                        ),
                        CodeBlock(str(sdata)),
                    ],
                )
            ],
        )

    def cell_section(self):
        shapes_key, gdf = get_boundaries(self.sdata, return_key=True)
        coord_system = get_intrinsic_cs(sdata, shapes_key)

        fig = plt.figure()
        sns.kdeplot(gdf.area)
        plt.xlabel("Area (coordinate_system_unit ^ 2)")

        return Section(
            "Cells",
            [
                SubSection(
                    "Number",
                    Paragraph(f"Number of cells:<br>{Message(f'{self.sdata.table.n_obs} cells')}"),
                ),
                SubSection(
                    "Areas",
                    [
                        Paragraph(
                            f"The cells areas are obtained based on the coordinate system '{coord_system}' for the '{shapes_key}' boundaries"
                        ),
                        Columns([Image(fig)]),
                    ],
                ),
            ],
        )

    def channel_section(self):
        image = get_spatial_image(self.sdata)

        subsections = [
            SubSection(
                "Names",
                Paragraph(
                    f"Channels names:<br>{Message(', '.join(list(image.coords['c'].values)))}"
                ),
            )
        ]

        if self._table_has(SopaKeys.UNS_HAS_INTENSITIES):
            fig = plt.figure()
            for channel, intensities in get_intensities(self.sdata).items():
                sns.kdeplot(intensities, label=channel)
            plt.xlabel("Intensity")
            plt.ylabel("Distribution density")

            subsections.append(
                SubSection(
                    "Per-cell intensity",
                    [
                        Paragraph(
                            "Each channel is averaged over each cell. We display the distribution of channel intensities over all cells:"
                        ),
                        Columns([Image(fig)]),
                    ],
                )
            )

        return Section("Channels", subsections)

    def transcripts_section(self):
        if not self._table_has(SopaKeys.UNS_HAS_TRANSCRIPTS):
            return None

        transcript_count = self.sdata.table.X.sum(0).A1
        low_average = transcript_count / self.sdata.table.n_obs < LOW_AVERAGE_COUNT

        QC_subsubsections = []
        if low_average.sum():
            QC_subsubsections.append(
                Paragraph(
                    f"{low_average.sum()} genes have a low count (less than {LOW_AVERAGE_COUNT} per cell, on average):"
                )
            )
            QC_subsubsections.append(
                Message(", ".join(self.sdata.table.var_names[low_average]), color="danger")
            )

        fig1 = plt.figure()
        sns.kdeplot(transcript_count)
        plt.xlabel("Total count per transcript")

        fig2 = plt.figure()
        sns.kdeplot(self.sdata.table.X.sum(1).A1)
        plt.xlabel("Transcript count per cell")

        QC_subsubsections.append(Columns([Image(fig1), Image(fig2)]))

        return Section("Transcripts", [SubSection("Quality controls", QC_subsubsections)])

    def representation_section(self):
        adata = self.sdata.table

        if self._table_has(SopaKeys.UNS_HAS_TRANSCRIPTS):
            sc.pp.log1p(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        colors = self._table_has(SopaKeys.UNS_CELL_TYPES, None)

        sc.pl.umap(adata, color=colors, show=False)
        fig = plt.gcf()

        return Section(
            "Representation",
            [
                SubSection("UMAP", Columns([Image(fig)])),
            ],
        )

    def compute_sections(self) -> list[Section]:
        sections = [
            self.general_section(),
            self.cell_section(),
            self.channel_section(),
            self.transcripts_section(),
            self.representation_section(),
        ]
        return [section for section in sections if section is not None]


sdata = spatialdata.read_zarr("/Volumes/T7_Quentin/data/test_sopa/dummy.zarr")
# sdata = None
write_report("sopa/io/report/test.html", sdata)
