import numpy as np

class FileHandlerXYZ:

    def __init__(self, file_name):
        self.file_name = file_name

    def read(self):

        symbols, x, y, z = np.loadtxt(self.file_name,
                                      skiprows=2,
                                      dtype={"names": ("atom", "x", "y", "z"),
                                             "formats": ("S2", "f4", "f4", "f4")},
                                      converters={3: _remove_zero_width_whitespace},
                                      encoding="UTF-8",
                                      unpack=True)

        symbols = symbols.astype(str)
        xyz = np.column_stack([x, y, z])

        return symbols, xyz

    def write(self, symbols, xyz, comment=""):

        np.savetxt(self.file_name,
                   np.c_[symbols, xyz],
                   encoding='utf-8',
                   fmt='%s %.10s %.10s %.10s',
                   header='\n'.join([str(xyz.shape[0]), comment]),
                   comments='')


def _remove_zero_width_whitespace(string):
    return string.strip(u"\u200b")



