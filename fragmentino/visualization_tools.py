#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
import numpy as np
import plotly.graph_objects as go

from fragmentino import Molecule


class Figure:
    def __init__(self, data):
        self.fig = go.Figure(data)
        self._update_layout()

    def _update_layout(self):
        self.fig.update_layout(
            showlegend=False,
            scene=dict(
                xaxis_title="",
                yaxis_title="",
                zaxis_title="",
                xaxis=dict(
                    showbackground=False,
                    tickvals=[],
                ),
                yaxis=dict(
                    showbackground=False,
                    tickvals=[],
                ),
                zaxis=dict(
                    showbackground=False,
                    tickvals=[],
                ),
            ),
            margin=dict(l=0, r=0, t=0, b=0),
        )

    def show_figure(self, **kwargs):
        self.fig.show(**kwargs)  # pragma: no cover

    def get_figure(self):
        return self.fig


class MoleculePlotter:
    def __init__(self, molecule, color=None):
        self.molecule = molecule
        self.color = color

    def get_bond_plot(self, label="bonds"):
        """Gets the data for plotting bonds using plotly

        Parameters
        ----------
        color : str
            String of color (hex or predefined color such as ``"pink"``).
            Default is ``color=None``, in which case the bonds are black.

        label : str
            Label assigned to the ploted atoms, default ``label="bonds"``

        Returns
        -------
        bonds : Scatter3d
            Data for 3D scatter plot (plotly)
        """
        if self.color == None:
            bond_color = "black"
        else:
            bond_color = self.color

        bonds = self.molecule.get_bonds()
        x_lines = []
        y_lines = []
        z_lines = []

        for bond in bonds:
            x_lines.append(self.molecule.xyz[bond[0], 0])
            x_lines.append(self.molecule.xyz[bond[1], 0])
            x_lines.append(None)

            y_lines.append(self.molecule.xyz[bond[0], 1])
            y_lines.append(self.molecule.xyz[bond[1], 1])
            y_lines.append(None)

            z_lines.append(self.molecule.xyz[bond[0], 2])
            z_lines.append(self.molecule.xyz[bond[1], 2])
            z_lines.append(None)

        bond_plot = go.Scatter3d(
            x=x_lines,
            y=y_lines,
            z=z_lines,
            name=label,
            mode="lines",
            line=dict(color=bond_color, width=10),
        )

        return bond_plot

    def get_atom_plot(self, label="atom"):
        """Gets the data for plotting atoms using plotly

        Parameters
        ----------
        color : str
            String of color (hex or predefined color such as "pink").
            Default is ``color=None``, in which case atoms are asigned color by their atomic number
            (https://en.wikipedia.org/wiki/CPK_coloring).

        label : str
            Label assigned to the ploted atoms, default ``label="atom"``

        Returns
        -------
        atoms : Scatter3d
            Data for 3D scatter plot (plotly)
        """
        from fragmentino.periodic_table import (
            Z_to_covalent_radius,
            Z_to_color,
        )

        marker_sizes = (
            20
            * np.fromiter(map(Z_to_covalent_radius, (self.molecule.Z)), dtype=float)
            * self.molecule.bond_factor
        )

        if self.color == None:  # Color by atomic number
            colors = []
            for Z in self.molecule.Z:
                colors.append(Z_to_color(Z))
        else:
            colors = self.color

        x = np.transpose(self.molecule.xyz[:, 0])
        y = np.transpose(self.molecule.xyz[:, 1])
        z = np.transpose(self.molecule.xyz[:, 2])

        atom_plot = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            name=label,
            marker=dict(
                size=marker_sizes,
                color=colors,
                opacity=1,
            ),
        )
        return atom_plot
