from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy.spatial import Delaunay
from shapely.geometry import LineString, Polygon


class Component:
    def __init__(self, adata: AnnData, delaunay: Delaunay, neighbors: np.ndarray) -> None:
        """Visitis one niche component to find border edges. Then, it assembles these border edges
        together to create the rings used as input to `shapely.geometry.Polygon`

        Args:
            adata: An `AnnData` object
            delaunay: The corresponding Delaunay graph
            neighbors: The neighbors of the delaunay graph
        """
        self.adata = adata
        self.delaunay = delaunay
        self.neighbors = neighbors

        self.lines: dict[int, list[int]] = {}  # values are ordered vertices indices
        self.limits: dict[int, tuple[int]] = {}  # values are (line_index, 0 if left / 1 if right)
        self._counter = 0

    def first_vertex_index(self) -> int:
        return next(iter(self.lines.values()))[0]

    def _standard_add(self, i: int, j: int, line_index_i: int, right_i: int):
        if right_i:
            del self.limits[self.lines[line_index_i][-1]]
            self.lines[line_index_i].append(j)
            self.limits[j] = (line_index_i, -1)
        else:
            del self.limits[i]
            self.lines[line_index_i].insert(0, j)
            self.limits[j] = (line_index_i, 0)

    def _insert_new(self, i: int, j: int):
        self._counter += 1
        self.lines[self._counter] = [i, j]
        self.limits[i] = (self._counter, 0)
        self.limits[j] = (self._counter, -1)

    def add_edge(self, i: int, j: int):
        line_index_i, right_i = self.limits.get(i, (None, None))
        line_index_j, right_j = self.limits.get(j, (None, None))

        if not line_index_i and not line_index_j:
            self._insert_new(i, j)

        elif line_index_i and not line_index_j:
            self._standard_add(i, j, line_index_i, right_i)

        elif line_index_j and not line_index_i:
            self._standard_add(j, i, line_index_j, right_j)

        else:  # join both lines
            del self.limits[j]
            del self.limits[i]

            if line_index_i == line_index_j:  # cycle
                self.lines[line_index_i].append(i if right_j else j)
                return

            line_i, line_j = self.lines[line_index_i], self.lines[line_index_j]
            del self.lines[line_index_j]
            self.limits[line_j[0 if right_j else -1]] = [line_index_i, right_i]

            if not (right_i ^ right_j):
                line_j = line_j[::-1]

            self.lines[line_index_i] = line_i + line_j if right_i else line_j + line_i

    def to_line(self, indices: list[int]) -> LineString:
        return LineString(self.adata.obsm["spatial"][indices])

    def __len__(self) -> int:
        return len(self.lines)

    @property
    def rings(self) -> list[LineString]:
        return [self.to_line(indices) for indices in self.lines.values()]

    @property
    def polygon(self) -> Polygon:
        rings = self.rings
        index_largest = np.argmax([Polygon(ring).area for ring in rings])

        return Polygon(
            rings[index_largest], rings[:index_largest] + rings[index_largest + 1 :]
        ).buffer(0)

    def visit(self, simplices_to_visit: set[int]):
        visited = set()
        component_simplices_to_visit = {simplices_to_visit.pop()}

        while component_simplices_to_visit:
            simplex_index = component_simplices_to_visit.pop()

            simplices_to_visit.discard(simplex_index)
            visited.add(simplex_index)

            simplex = self.delaunay.simplices[simplex_index]
            ngh_simplex_indices = self.neighbors[simplex_index]

            for i in range(3):
                ngh_simplex_index = ngh_simplex_indices[i]
                if ngh_simplex_index == -1:  # 'simplex' is at the boundary
                    self.add_edge(simplex[(i + 1) % 3], simplex[(i + 2) % 3])
                elif ngh_simplex_index not in visited:
                    component_simplices_to_visit.add(ngh_simplex_index)
