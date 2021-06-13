#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
from scipy.spatial import distance_matrix
import random


from fragmentino.molecule import Molecule
from fragmentino.periodic_table import Z_to_bond_length
from fragmentino import ContractableWeightedGraph
from fragmentino.visualization_tools import MoleculeFigure, MoleculePlotter


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
        self.g = ContractableWeightedGraph(max_fragment_size)
        self._fragment()

    def __getitem__(self, key):
        return self.g.vertices[key]

    def __iter__(self):
        for fragment in self.g.vertices:
            yield fragment

    def __repr__(self):
        return (
            f"{self.__class__.__name__}"
            + f" Fragments: {self.n_fragments}"
            + f" Capped bonds: {self.n_capped_bonds}"
        )

    @property
    def fragment_sizes(self):
        fragment_sizes = np.zeros(self.n_fragments, dtype=int)

        for i, fragment in enumerate(self):
            fragment_sizes[i] = fragment.size

        return fragment_sizes

    @property
    def n_fragments(self):
        return self.g.n_vertices

    @property
    def n_capped_bonds(self):
        return self.g.n_edges

    def write_separate(self, file_prefix):
        """Writes fragments to file. Fragment i is stored to ``file_prefix_fragment_i.xyz``

        Parameters
        ----------
        file_prefix : str
            prefix for file (with full or relative path).

        """
        for i, fragment in enumerate(self):
            fragment.write_xyz(file_prefix + "_fragment_" + str(i) + ".xyz")

    def write(self, file_prefix):
        """Writes fragments to a single file. Fragment i is stored to ``file_prefix_fragmented.xyz``

        Parameters
        ----------
        file_prefix : str
            prefix for file (with full or relative path).

        """
        m = Molecule(self.g.vertices[0].Z, self.g.vertices[0].xyz)
        for i, fragment in enumerate(self):
            if i > 0:
                m.merge(fragment)

        m.write_xyz(file_prefix + "_fragmented" + ".xyz")

    def add_H_to_capped_bonds(self):
        """
        Hydrogen is added with an apropriate bond length (given by covalent radii)
        in a direction given by the unit vector along the capped bond.
        """
        Z_H = 1
        for v1, v2 in self.g.edges:

            m1 = self.g.vertices[v1]
            m2 = self.g.vertices[v2]

            bonds = m1.get_bonds_to(m2)
            for bond in bonds:

                a1 = bond[0]
                a2 = bond[1]

                r = m2.xyz[a2, :] - m1.xyz[a1, :]
                n = r / (np.linalg.norm(r))

                # add H to m1
                length = Z_to_bond_length(m1.Z[a1], Z_H, m1.bond_factor)
                m1.add_atom(Z_H, m1.xyz[a1, :] + n * length)

                # add H to m2
                length = Z_to_bond_length(m2.Z[a2], Z_H, m2.bond_factor)
                m2.add_atom(Z_H, m2.xyz[a2, :] - n * length)

    def find_central_fragment(self):
        """
        Find central fragment by considering the center of
        mass of each fragment
        """
        CM = []
        for fragment in self:
            CM.append(fragment.center_of_mass)

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
        This process is :math:`\mathcal{O}(N^2)` scaling.
        """
        for i, fragment_i in enumerate(self):
            for j, fragment_j in enumerate(self):
                if j > i and fragment_i.same_size(fragment_j):
                    self.swap_fragments(i + 1, j)

    def plot_fragments(self, colors="random", **kwargs):
        """Plot fragments.

        Parameters
        ----------
        colors : str, optional
            Default is "random"
        kwargs
            Keyword arguments passed to :meth:`plotly.graph_objects.Figure.show`.
        """
        plots = []
        for fragment in self:

            if colors == "random":
                color = "#%06x" % random.randint(0, 0xFFFFFF)
            elif colors == "CPK":
                color = None

            plotter = MoleculePlotter(fragment, color)

            plots.append(plotter.get_bond_plot())
            plots.append(plotter.get_atom_plot())

        v = MoleculeFigure(data=plots)
        v.show(**kwargs)

    def _fragment(self):
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

        self.g.contract_by_smallest_weight()
