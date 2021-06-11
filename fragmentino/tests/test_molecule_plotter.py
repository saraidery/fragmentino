#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest
import os


from fragmentino import Molecule
from fragmentino import MoleculePlotter
from fragmentino import MoleculeFigure

class TestMoleculePlotter:
    def test_plotter(self):
        file_path = os.path.dirname(__file__)
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        plotter = MoleculePlotter(m)

        data=[plotter.get_atom_plot(), plotter.get_bond_plot()]
        v = MoleculeFigure(data=data)
        v.get_figure()

