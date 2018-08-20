#!/usr/bin/env python3

__author__ = "Pramit Barua"
__copyright__ = "Copyright 2018, INT, KIT"
__credits__ = ["Pramit Barua"]
__license__ = "INT, KIT"
__version__ = "1"
__maintainer__ = "Pramit Barua"
__email__ = ["pramit.barua@student.kit.edu", "pramit.barua@gmail.com"]


r'''
this method returns atom name and coordinate of atoms from xyz file
atom name is in list and coordinate is in numpy array
input argument:
    location: 
        directory of where the xyz file is
        example: C:\Users\Desktop
    file_name:
        the file name
        example: nt-4-0-3.xyz

output/return:
    atom name
    x, y and z coordinates
'''

import numpy as np
import os


def xyz_file_loader(location, file):
    status = True
    atoms = []
    coordinates = []

    file_name = os.path.join(location, file)

    if os.path.isfile(file_name):
        xyz = open(file_name)
        n_atoms = int(xyz.readline())
        title = xyz.readline()
        for line in xyz:
            try:
                atom, x, y, z = line.split()
                atoms.append(atom)
                coordinates.append([float(x), float(y), float(z)])
            except ValueError:
                pass
        xyz.close()

        if n_atoms != len(coordinates):
            status = False
            atoms = []
            coordinates = []
    else:
        status = False
        atoms = []
        coordinates = []

    return status, atoms, np.array(coordinates)


if __name__ == '__main__':
    status, atom, coordinate = xyz_file_loader(r'C:\Users\PRAMIT\Documents\MEGA\Master_thesis\pramit\python_code\related_to_thesis\zz_tube_xyz_generator', 'graphene_16_1.xyz')
    print(status)
