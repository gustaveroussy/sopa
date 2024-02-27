from __future__ import annotations

import base64
from io import BytesIO
from typing import Optional

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.figure import Figure

from .css import BULMA_CSS


class Renderable:
    """Object that can be transformed to string representing HTML"""

    @property
    def children(self) -> list["Renderable"]:
        if hasattr(self, "_children") and self._children is not None:
            if isinstance(self._children, list):
                return self._children
            return [self._children]
        return []

    @property
    def children_html(self) -> str:
        return "".join([str(child) for child in self.children])

    def children_rec(self) -> list["Renderable"]:
        return [self] + [cc for child in self.children for cc in child.children_rec()]

    def __str__(self) -> str:
        return self.children_html


class Title(Renderable):
    """Report title"""

    def __init__(self, text: str, level: int, subtitle: bool = False) -> None:
        self.text = text
        self.level = level
        self.subtitle = subtitle

    def __str__(self) -> str:
        return f"""<h1 class="{'subtitle' if self.subtitle else 'title'} is-{self.level}">{self.text}</h1>"""


class Paragraph(Renderable):
    """Report paragraph"""

    def __init__(self, text: str) -> None:
        self.text = text

    def __str__(self) -> str:
        return f"""<p>{self.text}</p>"""


class Message(Renderable):
    """Colored message"""

    def __init__(self, text: str, is_light: bool = True, color: str = "primary") -> None:
        self.text = text
        self.color = color
        self.is_light = is_light

    def __str__(self) -> str:
        return f"""
        <div class='notification is-{self.color} {'is-light' if self.is_light else ''}'>
            {self.text}
        </div>
        """


class Block(Renderable):
    """Block, i.e. padded div"""

    def __init__(self, content: list[Renderable]) -> None:
        self._children = content

    def __str__(self) -> str:
        return f"""<div class="block">{self.children_html}</div>"""


class CodeBlock(Renderable):
    """Block of code, like in the terminal"""

    def __init__(self, text: str) -> None:
        self.text = text

    def __str__(self) -> str:
        return f"""<pre>{self.text}</pre>"""


class ProgressBar(Renderable):
    """Progress bar"""

    def __init__(
        self,
        value: float,
        valuemax: int = 1,
        text: Optional[str] = None,
        color: str = "primary",
        is_light: bool = False,
    ) -> None:
        self.value = value
        self.valuemax = valuemax
        self.text = text
        self.color = color
        self.is_light = is_light

    def __str__(self) -> str:
        return f"""
        {"" if self.text is None else Paragraph(self.text)}
        <progress class="progress is-{self.color} {'is-light' if self.is_light else ''}"
            value="{self.value}"
            max="{self.valuemax}">{self.value}
        </progress>
    """


class Section(Renderable):
    """Section of the report"""

    def __init__(self, name: str, content: list["Section"] = None) -> None:
        self.name = name
        self._children = content
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
    """Sub-section of the report"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self) -> str:
        return f"""
        <section id="{self.id}" class="py-2">
            {Title(self.name, 4, subtitle=True)}
            {self.children_html}
        </section>
        """


class NavbarItem(Renderable):
    """One item in the nav bar"""

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
    """Left nav bar"""

    def __init__(self, sections: list[Section]) -> None:
        self._children = [NavbarItem(section) for section in sections]

    def __str__(self) -> str:
        return f"""
        {Title("Sopa report", 3)}
        {Message("This report was generated <br />by the <a href='https://github.com/gustaveroussy/sopa'>Sopa</a> HTML engine.")}
        {self.children_html}
    """


class Columns(Renderable):
    """Flex columns containers"""

    def __init__(self, content: list[Renderable], is_centered: bool = True) -> None:
        self._children = content
        self.is_centered = is_centered

    def __str__(self) -> str:
        return f"""
    <div block style="display: flex; {'justify-content: center;' if self.is_centered else ''}">
        {self.children_html}
    </div>
    """


class Image(Renderable):
    """Image renderer"""

    def __init__(
        self, fig: Figure, width: float = 50, extension: str = "png", pretty_legend: bool = True
    ):
        self.fig = fig
        self.width = width
        self.extension = extension
        self.pretty_legend = pretty_legend

    def make_figure_pretty(self):
        if self.pretty_legend and _has_handles(self.fig):
            self.fig.legend(
                bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0, frameon=False
            )
        sns.despine(fig=self.fig, offset=10, trim=True)

    def encod(self):
        self.make_figure_pretty()
        tmpfile = BytesIO()
        self.fig.savefig(tmpfile, format=self.extension, transparent=True, bbox_inches="tight")
        plt.close()
        return base64.b64encode(tmpfile.getvalue()).decode("utf-8")

    def __str__(self) -> str:
        return f"""<img src=\'data:image/{self.extension};base64,{self.encod()}\'  width="{self.width}%" height="auto"/>"""


def _has_handles(fig: Figure) -> bool:
    return any(len(ax.get_legend_handles_labels()[0]) for ax in fig.get_axes())


class Root(Renderable):
    """Whole report generator"""

    def __init__(self, sections: list[Section], doc_title: str = "Sopa report"):
        self.doc_title = doc_title
        self._children = sections
        self.nav = Navbar(sections)

    def sanity_check(self):
        section_ids = [section.id for section in self.children]
        assert len(section_ids) == len(set(section_ids)), "Sections IDs must be unique"

        subsections_ids = [sub.id for section in self.children for sub in section.children]
        assert len(subsections_ids) == len(set(subsections_ids)), "Subsections IDs must be unique"

    def write(self, path: str) -> None:
        self.sanity_check()

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
        <style>
        {BULMA_CSS}
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
