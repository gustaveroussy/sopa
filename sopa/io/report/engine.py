from pathlib import Path
from typing import Optional

import spatialdata
from spatialdata import SpatialData


class Renderable:
    @property
    def children(self) -> list["Renderable"]:
        if hasattr(self, "_children") and self._children is not None:
            return self._children
        return []

    @property
    def children_html(self) -> str:
        return "".join([str(child) for child in self.children])

    def children_rec(self) -> list["Renderable"]:
        return [self] + [cc for child in self.children for cc in child.children_rec()]

    def __str__(self) -> str:
        return self.children_html


class Paragraph(Renderable):
    def __init__(self, text: str) -> None:
        self.text = text

    def __str__(self) -> str:
        return f"""<p>{self.text}</p>"""


class CodeBlock(Renderable):
    def __init__(self, text: str) -> None:
        self.text = text

    def __str__(self) -> str:
        return f"""<pre>{self.text}</pre>"""


# class ProgressBar(Renderable):
#     def __init__(
#         self,
#         value: float,
#         th: list[float],
#         text: Optional[str] = None,
#         valuemax: int = 1,
#     ) -> None:
#         super().__init__(value, th)
#         self.text = text
#         self.valuemax = valuemax

#     @property
#     def perc(self):
#         return int(100 * self.value / self.valuemax)

#     def __str__(self) -> str:
#         return f"""
#     {"" if self.text is None else Paragraph(self.text)}
#     <div class="progress my-2">
#         <div
#             class="progress-bar bg-{self.color}"
#             role="progressbar"
#             style="width: {self.perc}%"
#             aria-valuenow="{self.value}"
#             aria-valuemin="0"
#             aria-valuemax="{self.valuemax}"
#           >
#             {self.perc}%
#         </div>
#     </div>
#     """


# class Alert(Renderable):
#     def __init__(self, value: float, th: list[float], text: str) -> None:
#         super().__init__(value, th)
#         self.text = text

#     def __str__(self) -> str:
#         return f"""
#     <div class="alert alert-{self.color}" role="alert">
#         {self.text}
#     </div>
#     """


class Section(Renderable):
    def __init__(self, name: str, children: list["Section"] = None) -> None:
        self.name = name
        self._children = children
        self.subtitle = False

    @property
    def id(self):
        return self.name.lower().replace(" ", "-")

    def __str__(self) -> str:
        return f"""
        <article class="message is-dark" id="{self.id}">
            <div class="message-header">
                <p>{self.name}</p>
            </div>
            <div class="message-body">
                {self.children_html}
            </div>
        </article>
        """


class SubSection(Section):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self) -> str:
        return f"""
        <section id="{self.id}" class="py-2">
            <h2 class="subtitle is-4">{self.name}</h2>
            {self.children_html}
        </section>
        """


class NavbarItem(Renderable):
    def __init__(self, section: Section) -> None:
        self.section = section

    def subsections(self):
        li = [
            f"<li><a href='#{subsection.id}'>{subsection.name}</a></li>"
            for subsection in self.section.children
        ]
        return "".join(li)

    def __str__(self) -> str:
        return f"""
    <p class="menu-label">{self.section.name}</p>
    <ul class="menu-list">
        {self.subsections()}
    </ul>
    """


class Navbar(Renderable):
    def __init__(self, sections: list[Section]) -> None:
        self._children = [NavbarItem(section) for section in sections]

    def __str__(self) -> str:
        return f"""
        <h3 class="title is-3">Sopa report</h3>
        <div class="notification is-primary is-light">
        This report was generated <br />by the
        <a href="https://github.com/gustaveroussy/sopa">Sopa</a> HTML
        engine.
        </div>
        {self.children_html}
    """


# class Flex(Renderable):
#     def __init__(self, children: list[Renderable]) -> None:
#         self._children = children

#     def __str__(self) -> str:
#         return f"""
#     <div class="d-flex p-2 justify-content-center">
#         {self.children}
#     </div>
#     """


# class Image(Renderable):
#     ASSETS_PATH: Path = None

#     def __init__(self, fig: Figure, name: str, width: float = 50, extension: str = "svg"):
#         self.fig = fig
#         self.name = name
#         self.width = width
#         self.extension = extension

#     @property
#     def fullname(self):
#         return f"{self.name}.{self.extension}"

#     def save(self):
#         self.ASSETS_PATH.mkdir(exist_ok=True, parents=True)
#         self.fig.savefig(self.ASSETS_PATH / self.fullname, bbox_inches="tight")

#     def __str__(self) -> str:
#         self.save()
#         return f"""
#     <img
#         src="{self.ASSETS_PATH.name}/{self.fullname}" width="{self.width}%" height="auto"
#     />"""


class Root(Renderable):
    CSS = "https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css"
    JS = "https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"
    ASSETS_NAME = "figures"

    def __init__(self, sections: list[Section], doc_title: str = "Sopa report"):
        self.doc_title = doc_title
        self._children = sections

        self.nav = Navbar(sections)

    def sanity_check(self):
        assert len(self.children) == len(
            set(section.id for section in self._children)
        ), "Sections IDs must be unique"

        # images = [child for child in self.children_rec() if isinstance(child, Image)]
        # assert len(images) == len(
        #     set(image.name for image in images)
        # ), "Images names must be unique"

    def write(self, path: str) -> None:
        self.sanity_check()

        # Image.ASSETS_PATH = Path(path).parent / self.ASSETS_NAME
        with open(path, "w") as f:
            f.write(str(self))

    def __str__(self) -> str:
        return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <title>{self.doc_title}</title>
        <link
        rel="stylesheet"
        href="https://cdn.jsdelivr.net/npm/bulma@0.9.4/css/bulma.min.css"
        />
        <style>
        .menu {{
            position: sticky;
            flex: 0 0 260px;
            overflow-y: auto;
            height: 100vh;
            top: 0;
        }}
        </style>
    </head>
    <body>
        <div class="is-flex is-flex-direction-row">
            <div class="mt-5 ml-5 menu">
                {self.nav}
            </div>
            <div class="p-5 block" style="flex: 1; overflow: hidden">
                {self.children_html}
            </div>
        </div>
    </body>
    </html>
    """


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
