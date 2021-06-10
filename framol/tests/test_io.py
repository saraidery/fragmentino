#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest
import os


from framol import io


class TestIO:
    def test_io_read(self):

        xyz_reference = np.array(
            [
                [0.86681, 0.60100, 0.00000],
                [-0.86681, 0.60144, 0.00000],
                [0.00000, -0.07579, 0.00000],
            ]
        )
        symbols_reference = ["H", "H", "O"]

        file_path = os.path.dirname(__file__)
        fh = io.FileHandlerXYZ(os.path.join(file_path, "small_molecule_1.xyz"))
        symbols, xyz = fh.read()

        assert np.allclose(xyz, xyz_reference)
        assert all(symbols == symbols_reference)

    def test_io_readwrite(self):

        symbols_reference = ["H", "He", "C"]
        xyz_reference = np.zeros((3, 3))

        file_path = os.path.dirname(__file__)
        fh = io.FileHandlerXYZ(os.path.join(file_path, "small_molecule_2.xyz"))
        fh.write(symbols_reference, xyz_reference)

        symbols, xyz = fh.read()

        assert np.allclose(xyz, xyz_reference)
        assert all(symbols == symbols_reference)
