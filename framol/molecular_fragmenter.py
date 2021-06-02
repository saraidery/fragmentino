import numpy as np
from scipy.spatial import distance_matrix


from framol.molecule import Molecule
from framol.periodic_table import covalent_radii
from framol import WeightedGraph


class MolecularFragmenter:
    """Handles the fragmentation of a molecule:

    - Makes a fragment for each atom. These atoms are the initial vertices of a graph
    - The bonds between atoms are edges for the graph
    - Contracts the graph by contracting over edges with smallest weights (shortest bonds)
    """

    def __init__(self, max_fragment_size, file_name):
        """Creates Molecular fragmenter

        Parameters
        ----------
        file_name : str
           Name xyz file to read (with full or relative path).
        """

        self.m = Molecule.from_xyz_file(file_name)
        self.g = WeightedGraph(max_fragment_size)

    def fragment(self):
        """Fragments the molecule"""
        for n, coord in zip(self.m.Z, self.m.xyz):
            m = Molecule(n, coord)
            self.g.add_vertex(m)

        bonds = self.m.get_bonds()
        for a1, a2, bond_length in bonds:
            self.g.add_edge(a1, a2, bond_length)

        self._merge_fragments()

    def store_fragments(self, file_prefix):
        """Writes fragments to file. Fragment i is stored to ``file_prefix_fragment_i.xyz``

        Parameters
        ----------
        file_prefix : str
            prefix for file (with full or relative path).

        """
        for i, molecule in enumerate(self.g.vertices):
            molecule.write(file_prefix + "_fragment_" + str(i) + ".xyz")

    def store_full(self, file_prefix):

        m = Molecule(self.g.vertices[0].Z, self.g.vertices[0].xyz)
        for i, molecule in enumerate(self.g.vertices):
            if i > 0:
                m.merge(molecule)

        m.write(file_prefix + "_full" + ".xyz")

    def add_H_to_capped_bonds(self):
        """
        Hydrogen is added with an apropriate bond length (given by covalent radii)
        in a direction given by the unit vector along the capped bond.
        """
        for v1, v2 in self.g.edges:

            m1 = self.g.vertices[v1]
            m2 = self.g.vertices[v2]

            bonds = m1.bonds_to(m2)

            for a1, a2, r in bonds:

                length = np.linalg.norm(r)
                n = r / length

                bond_to_H = covalent_radii[m1.Z[a1] - 1] + covalent_radii[0]
                m1.add_atom(1, m1.xyz[a1, :] + n * bond_to_H)

                bond_to_H = covalent_radii[m2.Z[a2] - 1] + covalent_radii[0]
                m2.add_atom(1, m2.xyz[a2, :] - n * bond_to_H)

    def _merge_fragments(self):
        """Merge fragments"""
        self.g.contract_by_smallest_weight()

    def find_center_fragment(self):
        """
        Find central fragment by considering the center of
        mass of each fragment
        """
        CM = []
        for vertex in self.g.vertices:
            CM.append(vertex.center_of_mass)

        CM = np.array(CM)
        i = np.linalg.norm(CM - np.mean(CM, axis=0), axis=1).argmin()
        return i

    def swap_fragments(self, v1, v2):
        """
        Swaps the order of two fragments
        """
        self.g.swap_vertices(v1, v2)
