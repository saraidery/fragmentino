#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
from scipy.spatial import distance_matrix


from fragmentino.io import FileHandlerXYZ
from fragmentino.periodic_table import (
    symbol_to_Z,
    Z_to_symbol,
    Z_to_covalent_radius,
    Z_to_atomic_weight,
    Z_to_color,
)


class Molecule:
    """Stores the molecule and its properties"""

    def __init__(self, Z, xyz, bond_factor=1.3):
        """Creates a molecule

        Parameters
        ----------
        Z : numpy.ndarray
            Atomic numbers (index n perodic table)

        xyz : numpy.ndarray
            Cartesian coordinates in Angstrom

        bond_factor : float
            Factor used to determine bonds. Default is ``bond_factor=1.3`` according to
            J. Chem. Phys. 117, 9160 (2002); https://doi.org/10.1063/1.1515483

        """
        self.xyz = np.atleast_2d(xyz)
        self.Z = np.atleast_1d(Z)
        self.bond_factor = bond_factor

    @classmethod
    def from_xyz_file(cls, file_name, bond_factor=1.3):
        """Creates a molecule by reading an xyz-file.

        Parameters
        ----------
        file_name : str
            File name of xyz-file with full or relative path.
        bond_factor : float
            Factor used to determine bonds. Default is ``bond_factor=1.3`` according to
            J. Chem. Phys. 117, 9160 (2002); https://doi.org/10.1063/1.1515483

        Returns
        -------
        molecule : Molecule

        """
        fh = FileHandlerXYZ(file_name)
        symbols, xyz = fh.read()
        Z = np.fromiter(map(symbol_to_Z, np.atleast_1d(symbols)), dtype=int)
        return cls(Z, xyz, bond_factor)

    @classmethod
    def from_molecules(cls, m1, m2, bond_factor=1.3):
        """Creates a molecule from two existing molecules

        Parameters
        ----------
        m1 : Molecule
            First molecule
        m2 : Molecule
            Second molecule
        bond_factor : float
            Factor used to determine bonds. Default is ``bond_factor=1.3`` according to
            J. Chem. Phys. 117, 9160 (2002); https://doi.org/10.1063/1.1515483

        Returns
        -------
        molecule : Molecule

        """

        Z = np.hstack((m1.Z, m2.Z))
        xyz = np.vstack((m1.xyz, m2.xyz))
        return cls(Z, xyz, bond_factor)

    def __repr__(self):
        return f"{self.__class__.__name__} {self.size}"

    def __len__(self):
        return self.size

    def __getitem__(self, key):
        return Molecule(self.Z[key], self.xyz[key])

    def __iter__(self):
        for Z, xyz in zip(self.Z, self.xyz):
            yield Z, xyz

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
        """Computes center of mass of molecule

        Returns
        -------
        CM : float
            Center of mass coordinate of molecule
        """
        atomic_weights = np.fromiter(map(Z_to_atomic_weight, (self.Z)), dtype=float)
        weighted_xyz = self.xyz * atomic_weights[:, None]

        CM = np.sum(weighted_xyz, axis=0)
        CM = CM / np.sum(atomic_weights)

        return CM

    @property
    def symbols(self):
        return list(map(Z_to_symbol, self.Z))

    def write_xyz(self, file_name, comment=""):
        """Writes the molecular geometry to an xyz-file.

        Parameters
        ----------
        file_name : str
            File name with full or relative path
        """
        fh = FileHandlerXYZ(file_name)
        fh.write(self.symbols, self.xyz, comment=comment)

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
        bond_length : ndarray
            Array storing sums of covalent radii.
            To be used to determine which atoms are bonded.

        """

        bond_length = np.zeros((self.size, self.size))
        for i, Z in enumerate(self.Z):
            bond_length[i, :] += Z_to_covalent_radius(Z)
            bond_length[:, i] += Z_to_covalent_radius(Z)

        bond_length = bond_length * self.bond_factor

        return bond_length

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
            List of bonds, given as ``[atom_1_index, atom_2_index, distance]``.
        """
        theoretical_bond_lengths = self._get_theoretical_covalent_bond_lengths()
        distances = self.distances

        rows, cols = np.where(distances < theoretical_bond_lengths)

        bonds = []
        for row, col in zip(rows, cols):
            if row < col:
                bonds.append([row, col, distances[row, col]])

        return bonds

    def get_bonds_to(self, other):
        """Determines bonds to another molecule.

        Parameters
        ----------
        other : Molecule
            The other molecule to merge with

        Returns
        -------
        bonds : list
            List of bonds, given as ``[atom_1_index, atom_2_index, distance]``.
        """
        distances = distance_matrix(self.xyz, other.xyz)

        theoretical_bond_lengths = np.zeros((self.size, other.size))
        for i, Z in enumerate(self.Z):
            theoretical_bond_lengths[i, :] += Z_to_covalent_radius(Z) * self.bond_factor

        for i, Z in enumerate(other.Z):
            theoretical_bond_lengths[:, i] += Z_to_covalent_radius(Z) * self.bond_factor

        rows, cols = np.where(distances < theoretical_bond_lengths)

        bonds = []
        for row, col in zip(rows, cols):
            bonds.append([row, col, distances[row, col]])

        return bonds

    def same_size(self, other):
        """Checks if two molecules are of the same size"""
        return self.size == other.size
