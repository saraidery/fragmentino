import numpy as np
import pytest


from framol import Molecule


class TestMolecule:

	@pytest.mark.parametrize("shape", [(3,3,4), (3,4)])
	def test_init_fail_coordinates(self, shape):

		xyz = np.zeros(shape)
		atomic_numbers = np.arange(1,4)

		with pytest.raises(ValueError, match="atomic coordinates"):

			m = Molecule(xyz, atomic_numbers)

	def test_init_fail_atomic_numbers(self):

		xyz = np.zeros((3,3))
		atomic_numbers = np.zeros((3,3))

		with pytest.raises(ValueError, match="atomic numbers"):

			m = Molecule(xyz, atomic_numbers)

	def test_init_fail_n_atoms(self):

		xyz = np.zeros((5,3))
		atomic_numbers = np.arange(1,4)

		with pytest.raises(ValueError, match="different number of atoms"):

			m = Molecule(xyz, atomic_numbers)

