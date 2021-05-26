import numpy as np


from framol.molecule import Molecule
from framol import WeightedGraph


class MolecularFragmenter:
    """Molecular fragmenter class"""

    def __init__(self, file_name, max_fragment_size):
        """Creates Molecular fragmenter"""
        self.m = Molecule.from_xyz_file(file_name)
        self.max_fragment_size = max_fragment_size

    def prepare_initial_fragments(self):

        fragments = []
        for n, coord in zip(self.m.atomic_numbers, self.m.xyz):

            m = Molecule(n, coord)
            fragments.append(m)

        self.g = WeightedGraph(fragments, self.max_fragment_size)

        bonds = self.m.get_covalent_bond_distances()
        distances = self.m.distances

        #is_bond = bonds > distances

        rows = np.where(distances < bonds)[0]
        cols = np.where(distances < bonds)[1]

        for row, col in zip(rows, cols):
            if (row < col):
                self.g.add_edge(row, col, distances[row, col])



