import numpy as np
import pytest
import os


from framol import MolecularFragmenter
from framol import Molecule


class TestFragmenter:
    def test_fragmentation(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(25, os.path.join(file_path, "small_molecule_1.xyz"))
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        f.fragment()

        assert np.allclose(np.sort(m.xyz, axis=0), np.sort(f.g.vertices[0].xyz, axis=0))
        assert np.allclose(np.sort(m.Z), np.sort(f.g.vertices[0].Z))

    def test_merge_fragments_2(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(30, os.path.join(file_path, "medium_molecule_1.xyz"))

        f.fragment()

        assert np.allclose([[0, 1]], f.g.edges[0])

    def test_merge_fragments_3(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "small_molecule_4.xyz"))

        f.fragment()

        edges = [[0, 2], [1, 2]]
        assert np.allclose(edges, f.g.edges)

    def test_add_H(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(2, os.path.join(file_path, "small_molecule_1.xyz"))

        f.fragment()

        f.add_H_to_capped_bonds()

        xyz_1 = [[-0.86681, 0.60144, 0.0000], [-0.23167467, 0.1052151, 0.0000]]

        xyz_2 = [
            [0.86681, 0.601, 0.0000],
            [0.0000, -0.07579, 0.0000],
            [-0.9936795, 0.7005619, 0.0000],
        ]

        assert np.allclose(xyz_1, f.g.vertices[0].xyz)
        assert np.allclose(xyz_2, f.g.vertices[1].xyz)

    def test_store(self):

        file_path = os.path.dirname(__file__)

        f = MolecularFragmenter(10, os.path.join(file_path, "small_molecule_1.xyz"))
        f.fragment()

        f.store_fragments(os.path.join(file_path, "small_molecule_1"))

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

        f.fragment()

        central_fragment_before = f.find_central_fragment()
        f.swap_fragments(central_fragment_before, 0)

        central_fragment_after = f.find_central_fragment()

        assert central_fragment_before != central_fragment_after
        assert central_fragment_after == 0

    def test_store_full(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(30, os.path.join(file_path, "medium_molecule_1.xyz"))
        f.fragment()

        f.store_full(os.path.join(file_path, "medium_molecule_1"))

        m1 = Molecule.from_xyz_file(
            os.path.join(file_path, "medium_molecule_1_full.xyz")
        )
        os.remove(os.path.join(file_path, "medium_molecule_1_full.xyz"))

        m2 = Molecule.from_xyz_file(os.path.join(file_path, "medium_molecule_1.xyz"))

        assert np.allclose(np.sort(m1.xyz, axis=0), np.sort(m2.xyz, axis=0))
        assert np.allclose(np.sort(m1.Z), np.sort(m2.Z))

    def test_group_fragments(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(10, os.path.join(file_path, "medium_molecule_1.xyz"))
        f.fragment()

        size, fragment_size_before = f.size
        f.group_fragments_by_size()
        fragment_size_after = f.size[1]

        before_reference = [ 4, 10,  9, 10]
        after_reference = [ 4, 10, 10,  9]
        assert np.allclose(fragment_size_before, before_reference)
        assert np.allclose(fragment_size_after, after_reference)

