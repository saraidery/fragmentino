import numpy as np
from scipy.spatial import distance_matrix


from framol.molecule import Molecule
from framol.periodic_table import covalent_radii
from framol import WeightedGraph


class MolecularFragmenter:
    """Molecular fragmenter class"""

    def __init__(self, file_name, max_fragment_size):
        """Creates Molecular fragmenter"""
        self.file_name = file_name
        self.m = Molecule.from_xyz_file(self.file_name)
        self.max_fragment_size = max_fragment_size

    def prepare_initial_fragments(self):
        """Prepare initial fragments"""
        fragments = []
        for n, coord in zip(self.m.atomic_numbers, self.m.xyz):

            m = Molecule(n, coord)
            fragments.append(m)

        self.g = WeightedGraph(fragments, self.max_fragment_size)

        bonds = self.m.get_covalent_bond_distances()
        distances = self.m.distances

        rows = np.where(distances < bonds)[0]
        cols = np.where(distances < bonds)[1]

        for row, col in zip(rows, cols):
            if row < col:
                self.g.add_edge(row, col, distances[row, col])

    def store_fragments(self):
        """Store fragments"""
        for i, molecule in enumerate(self.g.vertices):
            molecule.write(self.file_name[:-4] + "_fragment_" + str(i + 1) + ".xyz")

    def add_H_to_capped_bonds(self):
        """Add H to capped bonds"""

        # For each edge, search fragments for bonds between molecules
        # For each bond add H at position of other atom
        print(self.g.edges)
        for v1, v2 in self.g.edges:

            m1 = self.g.vertices[v1]
            m2 = self.g.vertices[v2]

            distances = distance_matrix(m1.xyz, m2.xyz)

            bonds = np.zeros((m1.size, m2.size))

            for i, Z in enumerate(m1.atomic_numbers):
                bonds[i, :] += covalent_radii[Z - 1]

            for i, Z in enumerate(m2.atomic_numbers):
                bonds[:, i] += covalent_radii[Z - 1]

            bonds = bonds * 1.1
            rows, cols = np.where(distances < bonds)

            for row, col in zip(rows, cols):

                # unit vector along capped bond
                r = (m2.xyz[col, :] - m1.xyz[row, :]) / distances[row, col]

                # Compute bond length by covalent radii
                bond_length_1 = 1.1 * (
                    covalent_radii[m1.atomic_numbers[row] - 1] + covalent_radii[0]
                )
                bond_length_2 = 1.1 * (
                    covalent_radii[m2.atomic_numbers[col] - 1] + covalent_radii[0]
                )

                m1.add_atom(1, m1.xyz[row, :] + r * bond_length_1)
                m2.add_atom(1, m2.xyz[col, :] - r * bond_length_2)
