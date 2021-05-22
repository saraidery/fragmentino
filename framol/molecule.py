class Molecule:
	"""Molecule class definition
	"""
	def __init__(self, coordinates, atomic_numbers):

		if coordinates.ndim != 2 and coordinates.shape[-1] != 3:
			raise ValueError('atomic coordinates must be of shape (N,3)')

		self.coordinates = coordinates
		self.atomic_numbers = atomic_numbers

	def __repr__(self):

		return f"{self.__class__.__name__} {self.size}"

	@property
	def size(self):
		return self.coordinates.shape[0]
