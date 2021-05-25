import numpy as np
import pytest
import os


from framol import Molecule
from framol.periodic_table import number_to_symbol, symbol_to_number
from framol import io


class TestMolecule:
    @pytest.mark.parametrize("shape", [(3, 3, 4), (3, 4)])
    def test_init_fail_xyz(self, shape):

        xyz = np.zeros(shape)
        atomic_numbers = np.arange(1, 4)

        with pytest.raises(ValueError, match="xyz must be"):

            m = Molecule(atomic_numbers, xyz)

    def test_init_fail_atomic_numbers(self):

        xyz = np.zeros((3, 3))
        atomic_numbers = np.zeros((3, 3))

        with pytest.raises(ValueError, match="atomic numbers"):

            m = Molecule(atomic_numbers, xyz)

    def test_init_fail_n_atoms(self):

        xyz = np.zeros((5, 3))
        atomic_numbers = np.arange(1, 4)

        with pytest.raises(ValueError, match="different number of atoms"):

            m = Molecule(atomic_numbers, xyz)

    def test_init(self):
        xyz = np.zeros((5, 3))
        atomic_numbers = np.arange(1, 6)

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
