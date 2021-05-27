import numpy as np
from scipy.spatial import distance_matrix


from framol.io import FileHandlerXYZ
from framol.periodic_table import number_to_symbol, symbol_to_number, covalent_radii


class Molecule:
    """Molecule class"""

    def __init__(self, atomic_numbers, xyz):
        """Creates a molecule

        Parameters
        ----------
        atomic_numbers : numpy array, int
            Atomic numbers (index n perodic table)

        xyz : numpy array, float
            Cartesian coordinates in Angstrom
        """
        self.xyz = np.array(xyz)
        self.atomic_numbers = np.array(atomic_numbers)

    @classmethod
    def from_xyz_file(cls, file_name):

        fh = FileHandlerXYZ(file_name)
        symbols, xyz = fh.read()
        atomic_numbers = np.fromiter(map(symbol_to_number, symbols), dtype=int)
        return cls(atomic_numbers, xyz)

    @classmethod
    def from_molecules(cls, m1, m2):

        atomic_numbers = np.hstack((m1.atomic_numbers, m2.atomic_numbers))
        xyz = np.vstack((m1.xyz, m2.xyz))

        return cls(atomic_numbers, xyz)

    def __repr__(self):
        return f"{self.__class__.__name__} {self.size}"

    @property
    def distances(self):
        """Distances between atoms in molecule"""
        return distance_matrix(self.xyz, self.xyz)

    @property
    def size(self):
        """Number of atoms in molecule"""
        return self.atomic_numbers.size

    def write(self, file_name):
        """Writes the molecular geometry to an xyz-file.

        Parameters
        ----------
        file_name : str
            File name with full or relative path
        """
        fh = FileHandlerXYZ(file_name)
        symbols = list(map(number_to_symbol, self.atomic_numbers))
        fh.write(symbols, self.xyz)

    def merge(self, other):
        """Merge

        Parameters
        ----------
        other : Molecule
            The other molecule to merge with

        """
        self.atomic_numbers = np.hstack(
            (self.atomic_numbers, other.atomic_numbers)
        )

        self.xyz = np.vstack((self.xyz, other.xyz))

    def get_covalent_bond_lengths(self):
        """Get covalent bond lengths

        Returns the matrix of covalent bond lengths given
        as the sum of the covalent radii of the atoms.

        Returns
        -------
        bonds : numpy.ndarray
            Array storing sums of covalent radii.
            To be used to determine which atoms are bonded.

        """

        bonds = np.zeros((self.size, self.size))
        for i, Z in enumerate(self.atomic_numbers):
            bonds[i, :] += covalent_radii[Z - 1]
            bonds[:, i] += covalent_radii[Z - 1]

        bonds = bonds
        return bonds

    def add_atom(self, atomic_number, xyz):
        """Add atom

        Appends an atom to the molecule

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
        self.atomic_numbers = np.append(self.atomic_numbers, atomic_number)
        self.xyz = np.vstack((self.xyz, xyz))
