#  fragmentino
#  Copyright (C) 2021 the authors of fragmentino

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np

# fmt: off
_periodic_table = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]
#https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
std_atomic_weight = [
    1.008, 4.003, 6.968, 9.012, 10.814, 12.011, 14.007, 15.999, 18.998, 20.1797,
    22.990, 24.3055, 26.982, 28.085, 30.974, 32.068, 35.452, 39.948, 39.098,
    40.078, 44.956, 47.867, 50.942,  51.996,  54.938, 55.845, 58.933, 58.693,
    3.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798,  85.468,
    87.62, 88.906, 91.224,  92.906, 95.95, 98.0, 101.07,  102.906, 106.42,
    107.868,  107.8682, 114.818, 118.710, 121.760, 127.60, 126.904,  131.293,
    132.905, 137.327, 138.905,  140.116, 140.907, 144.242, 145, 150.36, 151.964,
    157.25, 158.925, 162.500, 164.930, 167.259, 168.934,  173.054, 174.967,
    178.49, 180.948, 183.84, 186.207, 190.23, 192.217, 195.084, 196.967, 200.592, 204.384,  207.2,  208.980, 209, 210,
    222, 223, 226, 227, 232.038,  231.036, 238.029, 237, 244
]
covalent_radii = np.array([0.31, 0.28,
    1.28, 0.96, 0.84, 0.73, 0.71, 0.66, 0.57, 0.58,
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32,
    1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20,
    1.20, 1.16])

atom_color = ["#D2D2D2", "#00FFFF",   # H-He
              "#9933FF", "#009933", "#FF33CC", "#696969", "#0033FF", "#FF0000", "#00FF00", "#00FFFF",  # Li-Ne
              "#9933FF", "#009933", "#FF33CC", "#FF33CC", "#FF9900", "#FFFF00", "#00FF00", "#00FFFF",  # Na-Ar
              "#9933FF", "#009933", "#FF33CC", "#696969", "#FF33CC", "#FF33CC", "#FF33CC", "#FF6600",  # K-Fe
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#99000",  # Co-Br
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#6600CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC", "#FF33CC",
              "#FF33CC", "#FF33CC",
]
# fmt: on


def Z_to_symbol(Z):
    return _periodic_table[Z - 1]


def symbol_to_Z(symbol):
    n_elements = len(_periodic_table)
    periodic_table_dict = dict(
        zip(_periodic_table, np.arange(1, n_elements, dtype=int))
    )
    Z = periodic_table_dict[symbol]

    return Z


def Z_to_bond_length(Z1, Z2, scaling_factor):

    bond_length = (covalent_radii[Z1 - 1] + covalent_radii[Z2 - 1]) * scaling_factor

    return bond_length


def Z_to_covalent_radius(Z):

    return covalent_radii[Z - 1]


def Z_to_atomic_weight(Z):

    return std_atomic_weight[Z - 1]


def Z_to_color(Z):

    return atom_color[Z - 1]
