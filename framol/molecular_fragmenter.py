import numpy as np
from scipy.spatial import distance_matrix
import plotly.graph_objects as go
import random


from framol.molecule import Molecule
from framol.periodic_table import covalent_radii
from framol import WeightedGraph
from framol.visualization_tools import VisualizationTool


class MolecularFragmenter:
    """Handles the fragmentation of a molecule"""

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
        r"""Fragments the molecule in an :math:`\mathcal{O}(N^2)` procedure:

        - Makes a fragment for each atom. These atoms are the initial vertices of a graph
        - The bonds between atoms are edges for the graph
        - Contracts the graph by contracting over edges with smallest weights (shortest bonds)

        """
        for n, coord in zip(self.m.Z, self.m.xyz):
            m = Molecule(n, coord)
            self.g.add_vertex(m)

        bonds = self.m.get_bonds()
        for a1, a2, bond_length in bonds:
            self.g.add_edge(a1, a2, bond_length)

        self._merge_fragments()

    @property
    def size(self):
        return self.m.size

    @property
    def fragment_sizes(self):
        fragment_sizes = np.zeros(self.g.n_vertices, dtype=int)
        for i, fragment in enumerate(self.g.vertices):
            fragment_sizes[i] = fragment.size
        return fragment_sizes

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

                bond_to_H = (covalent_radii[m1.Z[a1] - 1] + covalent_radii[0])*m1.bond_factor
                m1.add_atom(1, m1.xyz[a1, :] + n * bond_to_H)

                bond_to_H = (covalent_radii[m2.Z[a2] - 1] + covalent_radii[0])*m2.bond_factor
                m2.add_atom(1, m2.xyz[a2, :] - n * bond_to_H)

    def _merge_fragments(self):
        """Merge fragments"""
        self.g.contract_by_smallest_weight()

    def find_central_fragment(self):
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

    def swap_fragments(self, f1, f2):
        """
        Swaps the order of two fragments

        Parameters
        ----------
        f1 : int
            Index of first fragment
        f2 : int
            Index of second fragment
        """
        self.g.swap_vertices(f1, f2)

    def group_fragments_by_size(self):
        r"""Groups fragments such that fragments of the same size follow each other

        Warning
        -------
        This should not be done before the molecule has been fragmented, since initially all fragments are of size 1
        and this process is :math:`\mathcal{O}(N^2)` scaling.
        """
        for i, vertex_i in enumerate(self.g.vertices):
            for j, vertex_j in enumerate(self.g.vertices):
                if j > i and vertex_i.same_size(vertex_j):
                    self.swap_fragments(i + 1, j)

    def plot_fragments(self, colors="by atom"):
        plots = []

        for i, molecule in enumerate(self.g.vertices):

            if (colors == "random"):
                color = "#%06x" % random.randint(i, 0xFFFFFF)
            elif (colors == "by atom"):
                color = None

            plots.append(molecule.get_atoms_plot(color))
            plots.append(molecule.get_bonds_plot(color))

        fig = go.Figure(data=plots)

        v = VisualizationTool()
        v.update_figure_layout(fig)
        v.show_figure(fig)
