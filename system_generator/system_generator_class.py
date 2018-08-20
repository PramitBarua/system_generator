#!/usr/bin/env python3

__author__ = "Pramit Barua"
__copyright__ = "Copyright 2018, INT, KIT"
__credits__ = ["Pramit Barua"]
__license__ = "INT, KIT"
__version__ = "1"
__maintainer__ = "Pramit Barua"
__email__ = ["pramit.barua@student.kit.edu", "pramit.barua@gmail.com"]

'''
task: 
check left_graphene_tube_overlap, graphene_unit_X, graphene_unit_Z is int

'''

import numpy as np
import os
import sys

sys.path.append(os.path.abspath('../NEGF_global'))

from xyz_file_loader import xyz_file_loader
from global_write import global_write


class SystemGenerator:
    def __init__(self, **kargs):
        self.a = 1.42
        # nearest neighbor distance

        system_block = []
        initial_coordinate = np.array([[0.0, 0.0, 0.0],
                                      [1.2297, 0.0, 0.71],
                                      [1.2297, 0.0, 2.13],
                                      [0.0, 0.0, 2.84]])
        status_X = True
        status_Z_left = True
        status_Z_right = True

        if 'input_parameter' in kargs:
            graphene_unit_X = kargs['input_parameter']['graphene']['num_unit_X']
            graphene_unit_Z = kargs['input_parameter']['graphene']['num_unit_Z']
            left_graphene_tube_overlap = kargs['input_parameter']['common']['left_graphene_tube_overlap']
            tube_unit = kargs['input_parameter']['tube']['tube_num_unit']
            right_graphene_tube_overlap = kargs['input_parameter']['common']['right_graphene_tube_overlap']

            half_graphene_unit_X = int(np.ceil(graphene_unit_X/2))

            graphene_dis_X = graphene_unit_X * np.sqrt(3) * self.a
            half_graphene_dis_X = half_graphene_unit_X * np.sqrt(3) * self.a
            # convert unit cell to angstrom

            tube_status, _, tube_coordinate = xyz_file_loader(
                kargs['input_parameter']['tube']['tube_dir'],
                kargs['input_parameter']['tube']['tube_name']
                )
            # read tube coordinates

            if tube_status:
                tube_radius = np.amax(tube_coordinate, axis=0)[1]
                status_X = self.check_distance(
                                graphene_dis_X,
                                2*tube_radius+10)
                # check the dimension of graphene sheet in X direction is
                # greater than or equal to the tube radius plus 10A
                # calculation is in angstrom

                status_Z_left = self.check_distance(
                                    graphene_unit_Z,
                                    left_graphene_tube_overlap)
                # checking, number of unit cell of graphene is greater than
                # or equal to the number of unit cell of tube that overlaps
                # with graphene sheet.
                # this check determines tube will not exceed graphene sheet.

                status_Z_right = self.check_distance(
                                    graphene_unit_Z,
                                    right_graphene_tube_overlap)
                # same concept as status_Z_left

                if status_X and status_Z_left and status_Z_right:
                    # shift tube in Y direction so that tube comes on top of
                    # the graphene sheet
                    tube_coordinate = self.add_value(
                                        tube_coordinate,
                                        tube_radius + kargs['input_parameter']['common']['graphene_tube_distance'],
                                        index=1)

                    # brings the tube in the middle of graphene sheet
                    tube_coordinate = self.add_value(
                                        tube_coordinate,
                                        half_graphene_dis_X,
                                        index=0)

                    # left graphene block
                    graphene_block = graphene_unit_Z - left_graphene_tube_overlap
                    for idx1 in range(graphene_block):
                        system_block = (system_block
                                        + self.graphene_sheet(
                                            initial_coordinate,
                                            0,
                                            graphene_unit_X))
                        initial_coordinate = self.add_value(
                                                initial_coordinate,
                                                3*self.a,
                                                2)

                        tube_coordinate = self.add_value(
                                            tube_coordinate,
                                            3*self.a,
                                            2)

                    # left graphene and tube overlap block
                    for idx1 in range(left_graphene_tube_overlap):
                        system_block = (system_block
                                        + self.graphene_sheet(
                                            initial_coordinate,
                                            0,
                                            half_graphene_unit_X))

                        system_block.append(tube_coordinate.tolist())

                        system_block = (system_block
                                        + self.graphene_sheet(
                                            initial_coordinate,
                                            half_graphene_unit_X+1,
                                            graphene_unit_X))

                        initial_coordinate = self.add_value(
                                                initial_coordinate,
                                                3*self.a,
                                                2)

                        tube_coordinate = self.add_value(
                                            tube_coordinate,
                                            3*self.a,
                                            2)

                    # middle tube block
                    for idx1 in range(tube_unit):
                        system_block.append(tube_coordinate.tolist())
                        initial_coordinate = self.add_value(
                                                initial_coordinate,
                                                3*self.a,
                                                2)
                        tube_coordinate = self.add_value(
                                            tube_coordinate,
                                            3*self.a,
                                            2)

                    # right graphene and tube overlap block
                    for idx1 in range(right_graphene_tube_overlap):
                        system_block = (system_block
                                        + self.graphene_sheet(
                                            initial_coordinate,
                                            0,
                                            half_graphene_unit_X))

                        system_block.append(tube_coordinate.tolist())

                        system_block = (system_block
                                        + self.graphene_sheet(
                                            initial_coordinate,
                                            half_graphene_unit_X,
                                            graphene_unit_X))
                        initial_coordinate = self.add_value(
                                                initial_coordinate,
                                                3*self.a,
                                                2)
                        tube_coordinate = self.add_value(
                                            tube_coordinate,
                                            3*self.a,
                                            2)

                    # right graphene block
                    if right_graphene_tube_overlap > 0:
                        # this 'if' defines that there is not overlap between
                        # graphene and tube in the right side
                        graphene_block = graphene_unit_Z - right_graphene_tube_overlap
                        for idx1 in range(graphene_block):
                            system_block = (system_block
                                            + self.graphene_sheet(
                                                initial_coordinate,
                                                0,
                                                graphene_unit_X))

                            initial_coordinate = self.add_value(
                                                    initial_coordinate,
                                                    3*self.a,
                                                    2)

                    system_block, cell_map = self.unit_cell(system_block)

                    file_name = 'coordinate ' + kargs['input_parameter']['common']['system_file_name']
                    self.write_file(kargs['input_parameter']['common']['output_file_dir'],
                                    file_name,
                                    system_block)

                    file_name = 'map ' + kargs['input_parameter']['common']['map_file_name']
                    self.write_file(kargs['input_parameter']['common']['output_file_dir'],
                                    file_name,
                                    cell_map)

            print(1)

    def write_file(self, location, file_name, array_element):
        file_data = [str(len(array_element))]
        file_data.append('Generated by Pramit Barua,'
                         'Master student in Karlsruhe Institute of Technology(KIT). '
                         'This code is written as a project of master thesis')
        for item in array_element:
            line_write = 'C '
            for idx, element in enumerate(item):
                if idx == len(item)-1:
                    line_write = line_write + str(element)
                else:
                    line_write = line_write + str(element) + ' '
            file_data.append(line_write)
        global_write(location, file_name, message = file_data)

    def unit_cell(self, system_block):
        count = 0
        map_block = []
        system = []
        for idx1 in range(len(system_block)):
            map_buf = []
            for idx2, item in enumerate(system_block[idx1]):
                system.append(item)
                count += 1
                map_buf.append(count)
            map_block.append(map_buf)
        return system, map_block

    def add_value(self, target, value, index):
        # add value in the ndarray columns.
        # target is a ndarray of numpy
        # value defines the amount to add (scalar)
        # index defines in which column the value will be added (int)
        target[:, index] += value
        return target

#     def elevate_Y(self, coordinates, ele_value):
#         coordinates[:, 1] += ele_value
#         return coordinates
# 
#     def shift_X(self, coordinates, shift_value):
#         coordinates[:, 0] += shift_value
#         return coordinates

    def check_distance(self, item_a, item_b):
        # check itemA >= itemB
        return True if item_a >= item_b else False

    def graphene_sheet(
            self, initial_coordinate, start_point, end_point, **kargs):
        # this method produce graphene sheet in X direction. Z direction it 
        # will produce always 1 unit

        end_point = end_point+1
        target = []
        for idx in range(start_point, end_point):
            distance = np.sqrt(3) * self.a * idx
            target.append(self.add_value(np.array(initial_coordinate), distance, 0).tolist())
        return target

#     def block_repeat_graphene(self):
#         pass
# 
#     def block_append(self):
#         pass
