import numpy as np
from scipy.spatial import distance_matrix


from framol.molecule import Molecule
from framol.periodic_table import covalent_radii
from framol import WeightedGraph


class MolecularFragmenter:
    """Molecular fragmenter class"""

    def __init__(self, max_fragment_size, file_name, path=""):
        """Creates Molecular fragmenter"""
        self.m = Molecule.from_xyz_file(path + file_name)
        self.max_fragment_size = max_fragment_size

    def prepare_initial_fragments(self):
        """Prepare initial fragments"""
        fragments = []
        for n, coord in zip(self.m.atomic_numbers, self.m.xyz):

            m = Molecule(n, coord)
            fragments.append(m)

        self.g = WeightedGraph(fragments, self.max_fragment_size)

        bonds = self.m.get_covalent_bond_lengths()
        distances = self.m.distances

        bonds = bonds

        rows = np.where(distances < bonds)[0]
        cols = np.where(distances < bonds)[1]

        for row, col in zip(rows, cols):
            if row < col:
                self.g.add_edge(row, col, distances[row, col])

    def store_fragments(self, file_prefix, path=""):
        """Store fragments"""
        for i, molecule in enumerate(self.g.vertices):
            molecule.write(path + file_prefix + "_fragment_" + str(i + 1) + ".xyz")

    def add_H_to_capped_bonds(self):
        """Add H to capped bonds"""

        for v1, v2 in self.g.edges:

            m1 = self.g.vertices[v1]
            m2 = self.g.vertices[v2]

            distances = distance_matrix(m1.xyz, m2.xyz)

            bonds = np.zeros((m1.size, m2.size))

            for i, Z in enumerate(m1.atomic_numbers):
                bonds[i, :] += covalent_radii[Z - 1]
            for i, Z in enumerate(m2.atomic_numbers):
                bonds[:, i] += covalent_radii[Z - 1]

            bonds = bonds
            rows, cols = np.where(distances < bonds)

            for row, col in zip(rows, cols):

                n = (m2.xyz[col, :] - m1.xyz[row, :]) / distances[row, col]

                bond_length_1 = (
                    covalent_radii[m1.atomic_numbers[row] - 1] + covalent_radii[0]
                )
                bond_length_2 = (
                    covalent_radii[m2.atomic_numbers[col] - 1] + covalent_radii[0]
                )

                m1.add_atom(1, m1.xyz[row, :] + n * bond_length_1)
                m2.add_atom(1, m2.xyz[col, :] - n * bond_length_2)

    def merge_fragments(self):
        """Merge fragments"""
        self.g.contract_by_smallest_weight()
