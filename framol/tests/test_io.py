import numpy as np
import pytest
import os


from framol import io


class TestMolecule:

   def test_io_read(self):

      xyz_reference =  np.array([[0.86681, 0.60144, 0.00000],
                                 [-0.86681, 0.60144, 0.00000],
                                 [0.00000 ,-0.07579, 0.00000]])
      symbols_reference = ['H', 'F', 'O']

      file_path = os.path.dirname(__file__)
      fh = io.FileHandlerXYZ(os.path.join(file_path, "small_molecule_1.xyz"))
      symbols, xyz = fh.read()

      assert np.allclose(xyz, xyz_reference)
      assert all(symbols == symbols_reference)

   def test_io_readwrite(self):

      symbols_reference = ['H', 'He', 'C']
      xyz_reference = np.zeros((3,3))

      file_path = os.path.dirname(__file__)
      fh = io.FileHandlerXYZ(os.path.join(file_path, "small_molecule_2.xyz"))
      fh.write(symbols_reference, xyz_reference)

      symbols, xyz = fh.read()

      assert np.allclose(xyz, xyz_reference)
      assert all(symbols == symbols_reference)
