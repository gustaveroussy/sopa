import numpy as np
from anndata import AnnData
from scipy.spatial import Delaunay
from shapely.geometry import LineString, Polygon


class Component:
    def __init__(self, adata: AnnData, delaunay: Delaunay, neighbors: np.array) -> None:
        self.adata = adata
        self.delaunay = delaunay
        self.neighbors = neighbors

        self.boundaries: dict[int, list[int]] = {}
        self.limits: dict[int, tuple[int]] = {}
        self._counter = 0

    def first_vertex_index(self) -> int:
        return next(iter(self.boundaries.values()))[0]

    def _standard_add(self, i: int, j: int, bound_index_i: int, loc_i: int):
        if loc_i == 0:
            del self.limits[i]
            self.boundaries[bound_index_i].insert(0, j)
            self.limits[j] = (bound_index_i, 0)
        else:
            del self.limits[self.boundaries[bound_index_i][-1]]
            self.boundaries[bound_index_i].append(j)
            self.limits[j] = (bound_index_i, -1)

    def _insert_new(self, i: int, j: int):
        self._counter += 1
        self.boundaries[self._counter] = [i, j]
        self.limits[i] = (self._counter, 0)
        self.limits[j] = (self._counter, -1)

    def add_edge(self, i: int, j: int):
        bound_index_i, loc_i = self.limits.get(i, (None, None))
        bound_index_j, loc_j = self.limits.get(j, (None, None))

        if bound_index_i and bound_index_j:
            del self.limits[j]
            del self.limits[i]

            if bound_index_i == bound_index_j:
                self.boundaries[bound_index_i].append(i if loc_j else j)
                return

            line_i = self.boundaries[bound_index_i]
            line_j = self.boundaries[bound_index_j]

            del self.boundaries[bound_index_j]

            if loc_i and not loc_j:  # ...i - j... > ...ij...
                self.limits[line_j[-1]] = [bound_index_i, -1]
                self.boundaries[bound_index_i] = line_i + line_j
                return

            if loc_j and not loc_i:  # i... - ...j > ...ji...
                self.limits[line_j[0]] = [bound_index_i, 0]
                self.boundaries[bound_index_i] = line_j + line_i
                return

            if loc_i and loc_j:  # ...i - ...j > ...ij...
                self.limits[line_j[0]] = [bound_index_i, -1]
                self.boundaries[bound_index_i] = line_i + line_j[::-1]
                return

            if not loc_i and not loc_j:  # i... - j... > ...ij...
                self.limits[line_i[-1]] = [bound_index_i, 0]
                self.limits[line_j[-1]] = [bound_index_i, -1]
                self.boundaries[bound_index_i] = line_i[::-1] + line_j
                return

        if bound_index_i and not bound_index_j:
            return self._standard_add(i, j, bound_index_i, loc_i)

        if bound_index_j and not bound_index_i:
            return self._standard_add(j, i, bound_index_j, loc_j)

        self._insert_new(i, j)

    def to_line(self, indices: list[int]) -> LineString:
        return LineString(self.adata.obsm["spatial"][indices])

    def __len__(self) -> int:
        return len(self.boundaries)

    @property
    def lines(self) -> list[LineString]:
        return [self.to_line(indices) for indices in self.boundaries.values()]

    @property
    def polygon(self) -> Polygon:
        lines = self.lines
        index_largest = np.argmax([Polygon(line).area for line in lines])

        return Polygon(
            lines[index_largest], lines[:index_largest] + lines[index_largest + 1 :]
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
