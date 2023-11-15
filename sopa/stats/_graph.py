import numpy as np
from anndata import AnnData
from shapely.geometry import LineString, Polygon


class Boundaries:
    def __init__(self, adata: AnnData) -> None:
        self.adata = adata

        self.boundaries: dict[int, list[int]] = {}
        self.limits: dict[int, tuple[int]] = {}
        self.n = 1

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
        self.boundaries[self.n] = [i, j]
        self.limits[i] = (self.n, 0)
        self.limits[j] = (self.n, -1)

        self.n += 1

    def add(self, i: int, j: int):
        (bound_index_i, loc_i), (bound_index_j, loc_j) = self.limits.get(
            i, (None, None)
        ), self.limits.get(j, (None, None))

        if bound_index_i and bound_index_j:
            del self.limits[j]
            del self.limits[i]

            if bound_index_i == bound_index_j:
                self.boundaries[bound_index_i].append(i if loc_j else j)
                # TODO: check limits not used anymore
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

    def __len__(self):
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
