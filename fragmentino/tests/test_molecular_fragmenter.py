#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest
import os


from fragmentino import MolecularFragmenter
from fragmentino import Molecule
from fragmentino import MoleculeFigure


class TestFragmenter:
    def test_fragmentation(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(25, os.path.join(file_path, "small_molecule_1.xyz"))
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        assert np.allclose(np.sort(m.xyz, axis=0), np.sort(f.g.vertices[0].xyz, axis=0))
        assert np.allclose(np.sort(m.Z), np.sort(f.g.vertices[0].Z))

    def test_merge_fragments_2(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(30, os.path.join(file_path, "medium_molecule_1.xyz"))

        assert np.allclose([[0, 1]], f.g.edges[0])

    def test_merge_fragments_3(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "small_molecule_4.xyz"))

        edges = [[0, 2], [1, 2]]
        assert np.allclose(edges, f.g.edges)

    def test_add_H(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(2, os.path.join(file_path, "small_molecule_1.xyz"))

        f.add_H_to_capped_bonds()

        xyz_1 = [[-0.86681, 0.60144, 0.0], [-0.29518825, 0.15483761, 0.0]]

        xyz_2 = [
            [0.86681, 0.601, 0.0],
            [0.0, -0.07579, 0.0],
            [-0.8943115, 0.62292665, 0.0],
        ]

        assert np.allclose(xyz_1, f.g.vertices[0].xyz)
        assert np.allclose(xyz_2, f.g.vertices[1].xyz)

    def test_write_separate(self):
        file_path = os.path.dirname(__file__)

        f = MolecularFragmenter(10, os.path.join(file_path, "small_molecule_1.xyz"))

        f.write_separate(os.path.join(file_path, "small_molecule_1"))

        m1 = Molecule.from_xyz_file(
            os.path.join(file_path, "small_molecule_1_fragment_0.xyz")
        )
        os.remove(os.path.join(file_path, "small_molecule_1_fragment_0.xyz"))

        m2 = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        assert np.allclose(np.sort(m1.xyz, axis=0), np.sort(m2.xyz, axis=0))
        assert np.allclose(np.sort(m1.Z), np.sort(m2.Z))

    def test_find_central_fragment_and_reorder(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(20, os.path.join(file_path, "medium_molecule_1.xyz"))

        central_fragment_before = f.find_central_fragment()
        f.swap_fragments(central_fragment_before, 0)

        central_fragment_after = f.find_central_fragment()

        assert central_fragment_before != central_fragment_after
        assert central_fragment_after == 0

    def test_write(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(30, os.path.join(file_path, "medium_molecule_1.xyz"))

        f.write(os.path.join(file_path, "medium_molecule_1"))

        m1 = Molecule.from_xyz_file(
            os.path.join(file_path, "medium_molecule_1_fragmented.xyz")
        )
        os.remove(os.path.join(file_path, "medium_molecule_1_fragmented.xyz"))

        m2 = Molecule.from_xyz_file(os.path.join(file_path, "medium_molecule_1.xyz"))

        print(np.sort(m1.xyz, axis=0))
        assert np.allclose(np.sort(m1.xyz, axis=0), np.sort(m2.xyz, axis=0))
        assert np.allclose(np.sort(m1.Z), np.sort(m2.Z))

    def test_group_fragments(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))

        sizes_before = f.fragment_sizes
        f.group_fragments_by_size()
        sizes_after = f.fragment_sizes

        before_reference = [4, 10, 9, 10]
        after_reference = [4, 10, 10, 9]
        assert np.allclose(sizes_before, before_reference)
        assert np.allclose(sizes_after, after_reference)

    def test_n_fragments(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))
        assert f.n_fragments == 4

    def test_n_capped_bonds(self):
        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))
        assert f.n_capped_bonds == 3

    def test_plot(self, monkeypatch):
        def mockreturn(v):
            return None

        monkeypatch.setattr(MoleculeFigure, "show", mockreturn)

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))
        f.plot_fragments()

    def test_plot_with_color(self, monkeypatch):
        def mockreturn(v):
            return None

        monkeypatch.setattr(MoleculeFigure, "show", mockreturn)

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))

        f.plot_fragments("CPK")

    def test_print(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))

        assert repr(f) == "MolecularFragmenter Fragments: 4 Capped bonds: 3"

    def test_indexing(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))

        assert repr(f[2]) == "Molecule 9"

    def test_order_fragments_by_centrality(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))

        f.order_fragments_by_centrality()

        CM = []
        for fragment in f:
            CM.append(fragment.center_of_mass)

        CM = np.array(CM)
        order = np.argsort(np.linalg.norm(CM - np.mean(CM, axis=0), axis=1))

        assert np.allclose(order, [0, 1, 2, 3])
