import numpy as np
import pytest
import os


from framol import Molecule
from framol.periodic_table import number_to_symbol, symbol_to_number
from framol import io

class TestMolecule:
    def test_init(self):
        xyz = np.zeros((5, 3))
        atomic_numbers = np.arange(1, 6)

        m = Molecule(atomic_numbers, xyz)

        assert np.allclose(xyz, m.xyz)
        assert np.allclose(atomic_numbers, m.atomic_numbers)

    def test_init_single_atom(self):
        xyz = np.array([0.0000, 0.0000, 5.0000])
        atomic_numbers = 2

        m = Molecule(atomic_numbers, xyz)

        assert np.allclose(xyz, m.xyz)
        assert np.allclose(atomic_numbers, m.atomic_numbers)

    def test_init_from_xyz_file(self):

        xyz_reference = np.array(
            [
                [0.86681, 0.60144, 0.00000],
                [-0.86681, 0.60144, 0.00000],
                [0.00000, -0.07579, 0.00000],
            ]
        )
        atomic_numbers_reference = [1, 9, 8]

        file_path = os.path.dirname(__file__)
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        assert np.allclose(xyz_reference, m.xyz)
        assert np.allclose(atomic_numbers_reference, m.atomic_numbers)

    def test_init_from_atoms(self):

        xyz_1 = [0.0000, 0.0000, 5.0000]
        atomic_numbers_1 = 2

        xyz_2 = np.array([0.0000, 5.0000, 0.0000])
        atomic_numbers_2 = 10

        m_1 = Molecule(atomic_numbers_1, xyz_1)
        m_2 = Molecule(atomic_numbers_2, xyz_2)

        m_3 = Molecule.from_molecules(m_1, m_2)

        xyz_reference = np.vstack((xyz_1, xyz_2))

        atomic_numbers_reference = np.array([2, 10])

        assert np.allclose(xyz_reference, m_3.xyz)
        assert np.allclose(atomic_numbers_reference, m_3.atomic_numbers)

    def test_init_from_molecules(self):

        xyz = np.zeros((2, 3))
        atomic_numbers = np.arange(1, 3)

        m1 = Molecule(atomic_numbers, xyz)

        xyz = np.ones((2, 3))
        atomic_numbers = np.arange(1, 3)

        m2 = Molecule(atomic_numbers, xyz)

        m3 = Molecule.from_molecules(m1, m2)

        xyz_reference =[[0.0000, 0.0000, 0.0000],
                        [0.0000, 0.0000, 0.0000],
                        [1.0000, 1.0000, 1.0000],
                        [1.0000, 1.0000, 1.0000]]
        atomic_numbers_reference = np.array([1, 2, 1, 2])


        assert np.allclose(xyz_reference, m3.xyz)
        assert np.allclose(atomic_numbers_reference, m3.atomic_numbers)

    def test_write(self):

        atomic_numbers_reference = np.array([1, 2, 6])
        xyz_reference = np.zeros((3, 3))

        m = Molecule(atomic_numbers_reference, xyz_reference)

        file_path = os.path.dirname(__file__)

        m.write(os.path.join(file_path, "small_molecule_2.xyz"))

        fh = io.FileHandlerXYZ(os.path.join(file_path, "small_molecule_2.xyz"))
        symbols, xyz = fh.read()

        atomic_numbers = np.fromiter(map(symbol_to_number, symbols), dtype=int)

        assert np.allclose(xyz_reference, xyz)
        assert np.allclose(atomic_numbers_reference, atomic_numbers)

    def test_size(self):
        xyz = np.zeros((5, 3))
        atomic_numbers = np.arange(1, 6)

        m = Molecule(atomic_numbers, xyz)
        size = m.size

        assert np.allclose(size, xyz.shape[0])

    def test_distances(self):

        file_path = os.path.dirname(__file__)
        m = Molecule.from_xyz_file(os.path.join(file_path, "small_molecule_1.xyz"))

        dstances_reference = np.array(
            [
                [0.0000000, 1.73362005, 1.10000005],
                [1.73362005, 0.0000000, 1.10000005],
                [1.10000005, 1.10000005, 0.0000000],
            ]
        )

        assert np.allclose(dstances_reference, m.distances)

    def test_merge_two_atoms(self):

        xyz_1 = [0.0000, 0.0000, 5.0000]
        atomic_numbers_1 = 2

        xyz_2 = np.array([0.0000, 5.0000, 0.0000])
        atomic_numbers_2 = 10

        m_1 = Molecule(atomic_numbers_1, xyz_1)
        m_2 = Molecule(atomic_numbers_2, xyz_2)

        m_1.merge(m_2)

        xyz_reference = np.vstack((xyz_1, xyz_2))

        atomic_numbers_reference = np.array([2, 10])

        assert np.allclose(xyz_reference, m_1.xyz)
        assert np.allclose(atomic_numbers_reference, m_1.atomic_numbers)

    def test_merge_two_molecules(self):

        xyz = np.zeros((2, 3))
        atomic_numbers = np.arange(1, 3)

        m1 = Molecule(atomic_numbers, xyz)

        xyz = np.ones((2, 3))
        atomic_numbers = np.arange(1, 3)

        m2 = Molecule(atomic_numbers, xyz)

        m1.merge(m2)

        xyz_reference =[[0.0000, 0.0000, 0.0000],
                        [0.0000, 0.0000, 0.0000],
                        [1.0000, 1.0000, 1.0000],
                        [1.0000, 1.0000, 1.0000]]
        atomic_numbers_reference = np.array([1, 2, 1, 2])


        assert np.allclose(xyz_reference, m1.xyz)
        assert np.allclose(atomic_numbers_reference, m1.atomic_numbers)

    def test_covalent_bond(self):
        xyz_1 = [0.0000, 0.0000, 5.0000]
        xyz_2 = [0.0000, 0.0000, 6.2000]
        xyz = np.vstack((xyz_1, xyz_2))
        atomic_numbers = np.array([7, 7])

        m = Molecule(atomic_numbers, xyz)

        bonds = m.get_covalent_bond_lengths()

        bonds_reference = [[1.42, 1.42],
                          [1.42, 1.42]]

        assert np.allclose(bonds_reference, bonds)

    def test_add_atom(self):
        xyz = np.zeros((2, 3))
        atomic_numbers = np.arange(1, 3)

        m = Molecule(atomic_numbers, xyz)

        xyz_atom = [0.0000, 0.0000, 5.0000]
        atomic_number_atom = 2

        xyz_reference =[[0.0000, 0.0000, 0.0000],
                        [0.0000, 0.0000, 0.0000],
                        [0.0000, 0.0000, 5.0000]]
        atomic_numbers_reference = np.append(atomic_numbers, atomic_number_atom)

        m.add_atom(atomic_number_atom, xyz_atom)

        assert np.allclose(xyz_reference, m.xyz)
        assert np.allclose(atomic_numbers_reference, m.atomic_numbers)

