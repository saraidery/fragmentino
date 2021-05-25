import numpy as np
from scipy.spatial import distance_matrix


from framol.io import FileHandlerXYZ
from framol.periodic_table import number_to_symbol, symbol_to_number


class Molecule:
    """Molecule class"""

    def __init__(self, atomic_numbers, xyz):
        """Creates a molecule

        Parameters
        ----------
        atomic_numbers : numpy array, int, shape(N)
            Atomic numbers (index n perodic table)

        xyz : numpy array, float, shape (N,3)
            Cartesian coordinates in Angstrom
        """
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ValueError("xyz must be of shape (N,3)")

        if atomic_numbers.ndim != 1:
            raise ValueError("atomic numbers must be of shape (N)")

        if atomic_numbers.shape[0] != xyz.shape[0]:
            raise ValueError("different number of atoms specified")

        self.xyz = xyz
        self.atomic_numbers = atomic_numbers

    @classmethod
    def from_xyz_file(cls, file_name):

        fh = FileHandlerXYZ(file_name)
        symbols, xyz = fh.read()
        atomic_numbers = np.fromiter(map(symbol_to_number, symbols), dtype=int)
        return cls(atomic_numbers, xyz)

    @classmethod
    def from_molecules(cls, m1, m2):

        atomic_numbers = np.concatenate((m1.atomic_numbers, m2.atomic_numbers))
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

        return self.xyz.shape[0]

    def write(self, file_name: str):
        """Writes the molecular geometry to an xyz-file.

        Parameters
        ----------
        file_name
            File name with full or relative path
        """
        fh = FileHandlerXYZ(file_name)
        symbols = list(map(number_to_symbol, self.atomic_numbers))
        fh.write(symbols, self.xyz)
