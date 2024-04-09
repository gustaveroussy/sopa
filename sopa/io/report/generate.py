from __future__ import annotations

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
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
    Root,
    Section,
    SubSection,
)

log = logging.getLogger(__name__)
warnings.simplefilter(action="ignore", category=FutureWarning)


def write_report(path: str, sdata: SpatialData):
    """Create a HTML report (or web report) after running Sopa.

    Note:
        This report is automatically generated based on a custom python-to-html engine

    Args:
        path: Path to the `.html` report that has to be created
        sdata: A `SpatialData` object, after running Sopa
    """
    sections = SectionBuilder(sdata).compute_sections()

    log.info(f"Writing report to {path}")
    Root(sections).write(path)


def _kdeplot_vmax_quantile(values: np.ndarray, quantile: float = 0.95):
    threshold = np.quantile(values, quantile)
    sns.kdeplot(values)
    plt.xlim(0, threshold)


class SectionBuilder:
    SECTION_NAMES = [
        "general_section",
        "cell_section",
        "channel_section",
        "transcripts_section",
        "representation_section",
    ]

    def __init__(self, sdata: SpatialData):
        self.sdata = sdata
        self.adata = self.sdata.tables.get(SopaKeys.TABLE)

    def _table_has(self, key, default=False):
        if SopaKeys.UNS_KEY not in self.adata.uns:
            return default
        return self.adata.uns[SopaKeys.UNS_KEY].get(key, default)

    def general_section(self):
        return Section(
            "General",
            [
                SubSection(
                    "SpatialData information",
                    [
                        Paragraph(
                            "Sopa is using <a href='https://spatialdata.scverse.org/en/latest/'>SpatialData</a> under the hood. This is how the object looks like:"
                        ),
                        CodeBlock(str(self.sdata)),
                    ],
                )
            ],
        )

    def cell_section(self):
        shapes_key, _ = get_boundaries(self.sdata, return_key=True)
        coord_system = get_intrinsic_cs(self.sdata, shapes_key)

        fig = plt.figure()
        _kdeplot_vmax_quantile(self.adata.obs[SopaKeys.AREA_OBS])
        plt.xlabel("Area (coordinate_system_unit ^ 2)")

        return Section(
            "Cells",
            [
                SubSection(
                    "Number",
                    Paragraph(f"Number of cells:<br>{Message(f'{self.adata.n_obs} cells')}"),
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
                    f"Channels names:<br>{Message(', '.join(map(str, list(image.coords['c'].values))))}"
                ),
            )
        ]

        if self._table_has(SopaKeys.UNS_HAS_INTENSITIES):
            fig = plt.figure()

            df_intensities = get_intensities(self.sdata)
            threshold = np.quantile(df_intensities.values.ravel(), 0.95)

            for channel, intensities in df_intensities.items():
                sns.kdeplot(intensities, label=channel)
            plt.xlim(0, threshold)
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

        mean_transcript_count = self.adata.X.mean(0).A1
        low_average = mean_transcript_count < LOW_AVERAGE_COUNT

        QC_subsubsections = []
        if low_average.sum():
            QC_subsubsections.append(
                Paragraph(
                    f"{low_average.sum()} genes have a low count (less than {LOW_AVERAGE_COUNT} per cell, on average):"
                )
            )
            QC_subsubsections.append(
                Message(", ".join(self.adata.var_names[low_average]), color="danger")
            )

        fig1 = plt.figure()
        _kdeplot_vmax_quantile(mean_transcript_count)
        plt.xlabel("Count per transcript (average / cells)")

        fig2 = plt.figure()
        _kdeplot_vmax_quantile(self.adata.X.sum(1).A1)
        plt.xlabel("Transcript count per cell")

        QC_subsubsections.append(Columns([Image(fig1), Image(fig2)]))

        return Section("Transcripts", [SubSection("Quality controls", QC_subsubsections)])

    def representation_section(self, max_obs: int = 400_000):
        if self._table_has(SopaKeys.UNS_HAS_TRANSCRIPTS):
            sc.pp.normalize_total(self.adata)
            sc.pp.log1p(self.adata)

        if self.adata.n_obs > max_obs:
            sc.pp.subsample(self.adata, n_obs=max_obs)

        log.info(f"Computing UMAP on {self.adata.n_obs} cells")

        sc.pp.pca(self.adata)
        sc.pp.neighbors(self.adata)
        sc.tl.umap(self.adata)

        colors = self._table_has(SopaKeys.UNS_CELL_TYPES, None)
        sc.pl.umap(self.adata, color=colors, show=False)

        return Section(
            "Representation",
            [
                SubSection("UMAP", Columns([Image(plt.gcf(), pretty_legend=False)])),
            ],
        )

    def compute_sections(self) -> list[Section]:
        sections = []

        for name in self.SECTION_NAMES:
            try:
                log.info(f"Writing {name}")
                section = getattr(self, name)()
                sections.append(section)
            except Exception as e:
                log.warn(f"Section {name} failed with error {e}")

        return [section for section in sections if section is not None]
