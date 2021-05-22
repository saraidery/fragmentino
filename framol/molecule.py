from framol.io import FileHandlerXYZ


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
		fh.write(self)

	#@classmethod
	#def from_xyz_file(cls, file_name):

		#fh = FileHandlerXYZ(file_name)
		#xyz, atomic_symbols = fh.read()
		# Atomic symbols -> atomic numbers
		#return cls(xyz, atomic_numbers)


