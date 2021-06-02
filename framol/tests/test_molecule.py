import numpy as np
import pytest
import os


from framol import Molecule
from framol.periodic_table import number_to_symbol, symbol_to_number
from framol import io


class TestMolecule:
    def test_init(self):
        xyz = np.zeros((5, 3))
        Z = np.arange(1, 6)

        m = Molecule(Z, xyz)

        assert np.allclose(xyz, m.xyz)
        assert np.allclose(Z, m.Z)

    def test_init_single_atom(self):
        xyz = np.array([0.0000, 0.0000, 5.0000])
        Z = 2

        m = Molecule(Z, xyz)

        assert np.allclose(xyz, m.xyz)
        assert np.allclose(Z, m.Z)

    def test_init_from_xyz_file(self):

        xyz_reference = np.array(
            [
                [0.86681, 0.60100, 0.00000],
                [-0.86681, 0.60144, 0.00000],
                [0.00000, -0.07579, 0.00000],
            ]
        )
        Z_reference = [1, 1, 8]

        file_path = os.path.dirname(__file__)
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        assert np.allclose(xyz_reference, m.xyz)
        assert np.allclose(Z_reference, m.Z)

    def test_init_from_atoms(self):

        xyz_1 = [0.0000, 0.0000, 5.0000]
        Z_1 = 2

        xyz_2 = np.array([0.0000, 5.0000, 0.0000])
        Z_2 = 10

        m_1 = Molecule(Z_1, xyz_1)
        m_2 = Molecule(Z_2, xyz_2)

        m_3 = Molecule.from_molecules(m_1, m_2)

        xyz_reference = np.vstack((xyz_1, xyz_2))

        Z_reference = np.array([2, 10])

        assert np.allclose(xyz_reference, m_3.xyz)
        assert np.allclose(Z_reference, m_3.Z)

    def test_init_from_molecules(self):

        xyz = np.zeros((2, 3))
        Z = np.arange(1, 3)

        m1 = Molecule(Z, xyz)

        xyz = np.ones((2, 3))
        Z = np.arange(1, 3)

        m2 = Molecule(Z, xyz)

        m3 = Molecule.from_molecules(m1, m2)

        xyz_reference = [
            [0.0000, 0.0000, 0.0000],
            [0.0000, 0.0000, 0.0000],
            [1.0000, 1.0000, 1.0000],
            [1.0000, 1.0000, 1.0000],
        ]
        Z_reference = np.array([1, 2, 1, 2])

        assert np.allclose(xyz_reference, m3.xyz)
        assert np.allclose(Z_reference, m3.Z)

    def test_write(self):

        Z_reference = np.array([1, 2, 6])
        xyz_reference = np.zeros((3, 3))

        m = Molecule(Z_reference, xyz_reference)

        file_path = os.path.dirname(__file__)

        m.write(os.path.join(file_path, "small_molecule_2.xyz"))

        fh = io.FileHandlerXYZ(os.path.join(file_path, "small_molecule_2.xyz"))
        symbols, xyz = fh.read()

        Z = np.fromiter(map(symbol_to_number, symbols), dtype=int)

        assert np.allclose(xyz_reference, xyz)
        assert np.allclose(Z_reference, Z)
        os.remove(os.path.join(file_path, "small_molecule_2.xyz"))

    def test_size(self):
        xyz = np.zeros((5, 3))
        Z = np.arange(1, 6)

        m = Molecule(Z, xyz)
        size = m.size

        assert np.allclose(size, xyz.shape[0])

    def test_distances(self):

        file_path = os.path.dirname(__file__)
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        dstances_reference = np.array(
            [
                [0.0000000, 1.7336201, 1.09972921],
                [1.7336201, 0.0000000, 1.10000005],
                [1.09972921, 1.10000005, 0.0000000],
            ]
        )

        assert np.allclose(dstances_reference, m.distances)

    def test_merge_two_atoms(self):

        xyz_1 = [0.0000, 0.0000, 5.0000]
        Z_1 = 2

        xyz_2 = np.array([0.0000, 5.0000, 0.0000])
        Z_2 = 10

        m_1 = Molecule(Z_1, xyz_1)
        m_2 = Molecule(Z_2, xyz_2)

        m_1.merge(m_2)

        xyz_reference = np.vstack((xyz_1, xyz_2))

        Z_reference = np.array([2, 10])

        assert np.allclose(xyz_reference, m_1.xyz)
        assert np.allclose(Z_reference, m_1.Z)

    def test_merge_two_molecules(self):

        xyz = np.zeros((2, 3))
        Z = np.arange(1, 3)

        m1 = Molecule(Z, xyz)

        xyz = np.ones((2, 3))
        Z = np.arange(1, 3)

        m2 = Molecule(Z, xyz)

        m1.merge(m2)

        xyz_reference = [
            [0.0000, 0.0000, 0.0000],
            [0.0000, 0.0000, 0.0000],
            [1.0000, 1.0000, 1.0000],
            [1.0000, 1.0000, 1.0000],
        ]
        Z_reference = np.array([1, 2, 1, 2])

        assert np.allclose(xyz_reference, m1.xyz)
        assert np.allclose(Z_reference, m1.Z)

    def test_add_atom(self):
        xyz = np.zeros((2, 3))
        Z = np.arange(1, 3)

        m = Molecule(Z, xyz)

        xyz_atom = [0.0000, 0.0000, 5.0000]
        atomic_number_atom = 2

        xyz_reference = [
            [0.0000, 0.0000, 0.0000],
            [0.0000, 0.0000, 0.0000],
            [0.0000, 0.0000, 5.0000],
        ]
        Z_reference = np.append(Z, atomic_number_atom)

        m.add_atom(atomic_number_atom, xyz_atom)

        assert np.allclose(xyz_reference, m.xyz)
        assert np.allclose(Z_reference, m.Z)

    def test_center_of_mass(self):
        """
        Computes the center of mass based on nuclear charges
        """

        file_path = os.path.dirname(__file__)
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        CM_reference = [0.00000, -2.80162898e-05, 0.0000]

        assert np.allclose(CM_reference, m.center_of_mass)

    def test_print(self):

        file_path = os.path.dirname(__file__)
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        assert repr(m) == "Molecule 3"
