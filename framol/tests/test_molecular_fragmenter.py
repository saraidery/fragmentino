import numpy as np
import pytest
import os


from framol import MolecularFragmenter
from framol import Molecule


class TestFragmenter:
    def test_prepare_initial_fragments(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(25, os.path.join(file_path, "small_molecule_1.xyz"))
        f.prepare_initial_fragments()

        edges = [[0, 2], [1, 2]]
        assert np.allclose(edges, f.g.edges)

    def test_merge_fragments(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(25, os.path.join(file_path, "small_molecule_1.xyz"))
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        f.prepare_initial_fragments()
        f.merge_fragments()

        assert np.allclose(np.sort(m.xyz, axis=0), np.sort(f.g.vertices[0].xyz, axis=0))
        assert np.allclose(np.sort(m.atomic_numbers), np.sort(f.g.vertices[0].atomic_numbers))


    def test_add_H(self):

        file_path = os.path.dirname(__file__)
        f = MolecularFragmenter(2, os.path.join(file_path, "small_molecule_1.xyz"))

        f.prepare_initial_fragments()
        f.merge_fragments()

        f.add_H_to_capped_bonds()

        xyz_1 = [[-0.86681, 0.60144, 0.0000],
                [-0.23167467, 0.1052151, 0.0000]]

        xyz_2 = [[ 0.86681, 0.601, 0.0000],
                 [ 0.0000, -0.07579, 0.0000],
                 [-0.9936795,  0.7005619, 0.0000]]

        assert np.allclose(xyz_1, f.g.vertices[0].xyz)
        assert np.allclose(xyz_2, f.g.vertices[1].xyz)
