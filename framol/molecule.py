import numpy as np


from framol.io import FileHandlerXYZ
from framol.periodic_table import number_to_symbol, symbol_to_number


class Molecule:
	"""Molecule class definition
	"""
	def __init__(self, coordinates, atomic_numbers):

		if coordinates.ndim != 2 or coordinates.shape[1] != 3:
			raise ValueError('atomic coordinates must be of shape (N,3)')

		if atomic_numbers.ndim != 1:
			raise ValueError('atomic numbers must be of shape (N)')

		if atomic_numbers.shape[0] != coordinates.shape[0]:
			raise ValueError('different number of atoms" specified')

		self.coordinates = coordinates
		self.atomic_numbers = atomic_numbers

	def __repr__(self):

		return f"{self.__class__.__name__} {self.size}"

	@property
	def size(self):

		return self.coordinates.shape[0]

	def write(self, file_name):

		fh = FileHandlerXYZ(file_name)
		symbols = list(map(number_to_symbol, self.atomic_numbers))
		fh.write(symbols, self.coordinates)

	@classmethod
	def from_xyz_file(cls, file_name):

		fh = FileHandlerXYZ(file_name)
		symbols, xyz = fh.read()
		atomic_numbers = np.fromiter(map(symbol_to_number, symbols), dtype=int)
		return cls(xyz, atomic_numbers)



