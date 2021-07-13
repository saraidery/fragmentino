#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import os


class FileHandlerXYZ:

    """Handles the reading and writing of xyz-files.

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
        self.file_name = os.path.expanduser(file_name.strip())

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
        with open(self.file_name, "w") as f:
            f.write(str(xyz.shape[0]) + "\n")
            f.write(comment + "\n")
            for symbol, pos in zip(symbols, xyz):
                f.write(
                    "{} {:15.10f} {:15.10f} {:15.10f}".format(
                        symbol, pos[0], pos[1], pos[2]
                    )
                )
                f.write("\n")


def _remove_zero_width_whitespace(string):
    """Remove non-zero whitespaces"""
    return string.strip(u"\u200b")
