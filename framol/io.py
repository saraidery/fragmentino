import numpy as np


class FileHandlerXYZ:

    """File handler xyz class
    Handles the readng and writing of xyz-files.

    The standard format of xyz-files assumed:
    line 1: Number of atoms
    line 2: Comment line
    line 3-: element  x-coordinate  y-coordinate  z-coordnate

    Units are Angstroms
    """

    def __init__(self, file_name):

        """Creates a xyz-file handler

        Parameters
        ----------
        file_name : str
            The name of the file, includng full or relative path.
        """

        self.file_name = file_name

    def read(self):

        """Read xyz-file"""

        symbols, x, y, z = np.loadtxt(
            self.file_name,
            skiprows=2,
            dtype={
                "names": ("atom", "x", "y", "z"),
                "formats": ("S2", "f4", "f4", "f4"),
            },
            converters={3: _remove_zero_width_whitespace},
            encoding="utf-8",
            unpack=True,
        )

        symbols = symbols.astype(str)
        xyz = np.column_stack([x, y, z])

        return symbols, xyz

    def write(self, symbols, xyz, comment=""):

        """Write xyz-file

        Parameters
        ----------
        symbols : list
            Symbols of the atoms in a molecule
        xyz : np.ndarray
            Cartesian coordinates of atoms of a molecule in Angstrom
        comment : str
            Optional comment for xyz-file comment line

        """

        np.savetxt(
            self.file_name,
            np.c_[symbols, xyz],
            encoding="utf-8",
            fmt="%s %.20s %.20s %.20s",
            header="\n".join([str(xyz.shape[0]), comment]),
            comments="",
        )


def _remove_zero_width_whitespace(string):

    """Remove non-zero whitespaces"""

    return string.strip(u"\u200b")
