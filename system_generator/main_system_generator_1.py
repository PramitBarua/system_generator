#!/usr/bin/env python3

__author__ = "Pramit Barua"
__copyright__ = "Copyright 2018, INT, KIT"
__credits__ = ["Pramit Barua"]
__license__ = "INT, KIT"
__version__ = "1"
__maintainer__ = "Pramit Barua"
__email__ = ["pramit.barua@student.kit.edu", "pramit.barua@gmail.com"]

'''

'''

import os
import sys
import argparse
import time
import numpy as np

sys.path.append(os.path.abspath('../NEGF_global'))

from yaml_file_loader import yaml_file_loader
from xyz_file_loader import xyz_file_loader
from global_write import global_write

from system_generator_class import SystemGenerator


if __name__ == '__main__':
    print('system generator start')
    status = True
    start_time = time.time()

    parser = argparse.ArgumentParser(description='Calculate and display DOS ' +
                                     'and transmission of the quantum system.')
    parser.add_argument("Folder_Name",
                        help="Name of the folder that contains " +
                        "'parameter_system_generator.yml' file")
    args = parser.parse_args()

    input_status, input_parameter = yaml_file_loader(args.Folder_Name,
                                                     'parameter_system_generator.yml')

    if input_status:
        system = SystemGenerator(input_parameter = input_parameter)
#         system.graphene_tube_left()
    else:
        print(' === Error: yaml file failed to load === ')
