import spatialdata
from spatialdata import SpatialData

from .engine import CodeBlock, Paragraph, Root, Section, SubSection


def write_report(path: str, sdata):
    section1 = Section(
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
    section2 = Section("Section 2!", [SubSection("subsection12"), SubSection("subsection22")])
    section3 = Section("Section 3!", [SubSection("subsection13"), SubSection("subsection23")])
    section4 = Section("Section 4!", [SubSection("subsection14"), SubSection("subsection24")])
    section5 = Section("Section 5!", [SubSection("subsection15"), SubSection("subsection25")])

    root = Root([section1, section2, section3, section4, section5])
    root.write(path)


# sdata = spatialdata.read_zarr("/Volumes/T7_Quentin/data/vizgen/test.zarr")
sdata = None
write_report("test.html", sdata)
