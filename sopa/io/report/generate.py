# import spatialdata
# from spatialdata import SpatialData
import matplotlib.pyplot as plt

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


def write_report(path: str, sdata):
    sections = SectionBuilder(sdata).compute_sections()

    Root(sections).write(path)


class SectionBuilder:
    def __init__(self, sdata):
        self.sdata = sdata

    def general_section(self):
        return Section(
            "General",
            [
                SubSection(
                    "SpatialData information",
                    [
                        Paragraph(f"Number of cells:<br>{Message('130 cells')}"),
                        Paragraph(
                            f"Sopa is using <a href='https://spatialdata.scverse.org/en/latest/'>SpatialData</a> under the hood. This is how the object looks like:"
                        ),
                        CodeBlock(str(sdata)),
                    ],
                )
            ],
        )

    def channel_section(self):
        return Section(
            "Channels",
            [
                SubSection("subsection121"),
                SubSection("subsection222", ProgressBar(0.4, text="Hello")),
            ],
        )

    def transcripts_section(self):
        return Section("Transcripts", [SubSection("subsection123"), SubSection("subsection224")])

    def annotation_section(self):
        fig = plt.figure()
        plt.scatter([1, 2], [3, 4])
        return Section(
            "Annotation",
            [
                SubSection("subsection125"),
                SubSection("subsection226", Columns([Image(fig), Image(fig)])),
            ],
        )

    def compute_sections(self) -> list[Section]:
        return [
            self.general_section(),
            self.channel_section(),
            self.transcripts_section(),
            self.annotation_section(),
        ]


# sdata = spatialdata.read_zarr("/Volumes/T7_Quentin/data/vizgen/test.zarr")
sdata = None
write_report("sopa/io/report/test.html", sdata)
