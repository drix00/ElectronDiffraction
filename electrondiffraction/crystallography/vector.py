#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: vector
   :synopsis: Vector operation in a crystal system.

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Vector operation in a crystal system.
"""

###############################################################################
# Copyright 2017 Hendrix Demers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###############################################################################

# Standard library modules.

# Third party modules.
import numpy as np

# Local modules.

# Project modules.

# Globals and constants variables.


def dot_product(crystal, vector_p, vector_q):
    g_ij_nm2 = crystal.gij_nm2

    vector1 = np.dot(g_ij_nm2, np.array(vector_q).transpose())
    magnitude = np.dot(np.array(vector_p), vector1)

    return magnitude


def distance(crystal, vector_p, vector_q):
    d = np.sqrt(dot_product(crystal, vector_p, vector_q))

    return d


def length(crystal, vector_p):
    return distance(crystal, vector_p, vector_p)


def distance_points(crystal, point_p, point_q):
    vector_d = point_p - point_q

    return distance(crystal, vector_d, vector_d)


def angle_rad(crystal, vector_p, vector_q):
    denominator = dot_product(crystal, vector_p, vector_q)
    norm_p = length(crystal, vector_p)
    norm_q = length(crystal, vector_q)

    factor = denominator / (norm_p * norm_q)

    angle_rad = np.arccos(factor)

    return angle_rad


def angle2_rad(crystal, vector_p, vector_q):
    matrix_left = np.vstack([vector_p, vector_q])
    matrix_right = matrix_left.transpose()

    matrix = np.dot(crystal.gij_nm2, matrix_right)
    matrix = np.dot(matrix_left, matrix)

    nominator = matrix[0, 1]
    denominator = np.sqrt(matrix[0, 0]) * np.sqrt(matrix[1, 1])
    factor = nominator / denominator

    angle_rad = np.arccos(factor)

    return angle_rad
