import numpy as np
from scipy.spatial import distance_matrix
import plotly.graph_objects as go


from framol.io import FileHandlerXYZ
from framol.periodic_table import (
    number_to_symbol,
    symbol_to_number,
    covalent_radii,
    std_atomic_weight,
    atom_color,
)


class Molecule:
    def __init__(self, Z, xyz, bond_factor=1.3):
        """Creates a molecule

        Parameters
        ----------
        Z : numpy array(int)
            Atomic numbers (index n perodic table)

        xyz : numpy array, float
            Cartesian coordinates in Angstrom
        """
        self.xyz = np.atleast_2d(xyz)
        self.Z = np.atleast_1d(Z)
        self.bond_factor=bond_factor

    @classmethod
    def from_xyz_file(cls, file_name):
        """Creates a molecule by reading an xyz-file.

        Parameters
        ----------
        file_name : str
            File name of xyz-file with full or relative path.
        """
        fh = FileHandlerXYZ(file_name)
        symbols, xyz = fh.read()
        Z = np.fromiter(map(symbol_to_number, np.atleast_1d(symbols)), dtype=int)
        return cls(Z, xyz)

    @classmethod
    def from_molecules(cls, m1, m2):
        """Creates a molecule from two existing molecules

        Parameters
        ----------
        m1 : Molecule
            First molecule
        m2 : Molecule
            Second molecule

        Returns
        -------
        molecule : Molecule
        """

        Z = np.hstack((m1.Z, m2.Z))
        xyz = np.vstack((m1.xyz, m2.xyz))
        return cls(Z, xyz)

    def __repr__(self):
        return f"{self.__class__.__name__} {self.size}"

    @property
    def distances(self):
        """Distances between atoms in molecule"""
        return distance_matrix(self.xyz, self.xyz)

    @property
    def size(self):
        """Number of atoms in molecule"""
        return self.Z.size

    @property
    def center_of_mass(self):

        atomic_weights = np.array(std_atomic_weight)[self.Z - 1]
        weighted_xyz = self.xyz * atomic_weights[:, None]
        CM = np.sum(weighted_xyz, axis=0)
        CM = CM / np.sum(atomic_weights)
        return CM

    def write(self, file_name):
        """Writes the molecular geometry to an xyz-file.

        Parameters
        ----------
        file_name : str
            File name with full or relative path
        """
        fh = FileHandlerXYZ(file_name)
        symbols = list(map(number_to_symbol, self.Z))
        fh.write(symbols, self.xyz)

    def merge(self, other):
        """Appends another molecule to it self.

        Parameters
        ----------
        other : Molecule
            The other molecule to merge with

        """
        self.Z = np.hstack((self.Z, other.Z))
        self.xyz = np.vstack((self.xyz, other.xyz))

    def _get_theoretical_covalent_bond_lengths(self):
        """Returns the matrix of sums of covalent radii for each pair
        of atoms in the molecule.

        Returns
        -------
        bonds : numpy.ndarray
            Array storing sums of covalent radii.
            To be used to determine which atoms are bonded.

        """

        bonds = np.zeros((self.size, self.size))
        for i, Z in enumerate(self.Z):
            bonds[i, :] += covalent_radii[Z - 1]*self.bond_factor
            bonds[:, i] += covalent_radii[Z - 1]*self.bond_factor

        bonds = bonds
        return bonds

    def add_atom(self, atomic_number, xyz):
        """Appends an atom to the molecule

        Parameters
        ----------
        atomic_number : int
            Atomic number of the appended atom
        xyz : numpy.ndarray
            Cartesian coordinates of appended atom

        Note
        ----

        Changes the instance of the molecule

        """
        self.Z = np.append(self.Z, atomic_number)
        self.xyz = np.vstack((self.xyz, xyz))

    def get_bonds(self):
        """Determines the bonds of the molecule

        Returns
        -------
        bonds : list
            List of bonds, given as ``[atom_idx_1, atom_idx_2, bond_length]``.
        """
        sum_covalent_radii = self._get_theoretical_covalent_bond_lengths()
        distances = self.distances

        rows, cols = np.where(distances < sum_covalent_radii)

        bonds = []
        for row, col in zip(rows, cols):
            if row < col:
                bonds.append([row, col, distances[row, col]])

        return bonds

    def bonds_to(self, other):
        """Determines bonds to another molecule.

        Parameters
        ----------
        other : Molecule
            The other molecule to merge with

        Returns
        -------
        bonds : list
            List of bonds, given as ``[atom_idx_1, atom_idx_2, r]``
            where r is the vector along the bond.
        """
        distances = distance_matrix(self.xyz, other.xyz)

        sum_covalent_radii = np.zeros((self.size, other.size))
        for i, Z in enumerate(self.Z):
            sum_covalent_radii[i, :] += covalent_radii[Z - 1]*self.bond_factor
        for i, Z in enumerate(other.Z):
            sum_covalent_radii[:, i] += covalent_radii[Z - 1]*self.bond_factor

        rows, cols = np.where(distances < sum_covalent_radii)

        bonds = []
        for row, col in zip(rows, cols):

            r = other.xyz[col, :] - self.xyz[row, :]
            bonds.append([row, col, r])

        return bonds

    def same_size(self, other):
        return self.size == other.size

    def plot(self, color=None):

        atom_bonds = self.get_bonds_plot()
        atoms = self.get_atoms_plot()

        fig = go.Figure(data=[atom_bonds, atoms])

        fig.update_layout(
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
        fig.show()

    def get_bonds_plot(self, color=None, label='bonds'):

        if (color == None):
            bond_color = "black"
        else:
            bond_color = color

        bonds = self.get_bonds()
        x_lines = []
        y_lines = []
        z_lines = []

        # create bonds
        for bond in bonds:
            x_lines.append(self.xyz[bond[0], 0])
            x_lines.append(self.xyz[bond[1], 0])
            x_lines.append(None)
            y_lines.append(self.xyz[bond[0], 1])
            y_lines.append(self.xyz[bond[1], 1])
            y_lines.append(None)
            z_lines.append(self.xyz[bond[0], 2])
            z_lines.append(self.xyz[bond[1], 2])
            z_lines.append(None)

        atom_bonds = go.Scatter3d(
            x=x_lines,
            y=y_lines,
            z=z_lines,
            name=label,
            mode="lines",
            line=dict(color=bond_color, width=10),
        )

        return atom_bonds


    def get_atoms_plot(self, color=None, label='atom'):

        sizes = 25 * covalent_radii[self.Z - 1]*self.bond_factor

        if (color == None):
            colors = []
            for Z in self.Z:
                colors.append(atom_color[Z - 1])
        else:
            colors = color

        x = np.transpose(self.xyz[:,0])
        y = np.transpose(self.xyz[:,1])
        z = np.transpose(self.xyz[:,2])

        atoms = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            name=label,
            marker=dict(
                size=sizes,
                color=colors,
                opacity=1,
            ),
        )

        return atoms
