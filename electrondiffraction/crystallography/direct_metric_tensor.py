#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: direct_metric_tensor
   :synopsis: Equations to compute the direct metric tensor for each crystal system.

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Equations to compute the direct metric tensor for each crystal system.
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


def gc_nm2(a_nm):
    tensor_nm2 = np.zeros((3, 3))

    tensor_nm2[0, 0] = a_nm * a_nm
    tensor_nm2[1, 1] = a_nm * a_nm
    tensor_nm2[2, 2] = a_nm * a_nm

    return tensor_nm2


def gt_nm2(a_nm, c_nm):
    tensor_nm2 = np.zeros((3, 3))

    tensor_nm2[0, 0] = a_nm * a_nm
    tensor_nm2[1, 1] = a_nm * a_nm
    tensor_nm2[2, 2] = c_nm * c_nm

    return tensor_nm2


def go_nm2(a_nm, b_nm, c_nm):
    tensor_nm2 = np.zeros((3, 3))

    tensor_nm2[0, 0] = a_nm * a_nm
    tensor_nm2[1, 1] = b_nm * b_nm
    tensor_nm2[2, 2] = c_nm * c_nm

    return tensor_nm2


def gh_nm2(a_nm, c_nm):
    tensor_nm2 = np.zeros((3, 3))

    tensor_nm2[0, 0] = a_nm * a_nm
    tensor_nm2[1, 1] = a_nm * a_nm
    tensor_nm2[2, 2] = c_nm * c_nm
    tensor_nm2[0, 1] = -a_nm * a_nm / 2.0
    tensor_nm2[1, 0] = -a_nm * a_nm / 2.0

    return tensor_nm2


def gr_nm2(a_nm, alpha_rad):
    tensor_nm2 = np.zeros((3, 3))
    cos_alpha = np.cos(alpha_rad)

    tensor_nm2[0, 0] = a_nm * a_nm
    tensor_nm2[0, 1] = a_nm * a_nm * cos_alpha
    tensor_nm2[0, 2] = a_nm * a_nm * cos_alpha
    tensor_nm2[1, 0] = a_nm * a_nm * cos_alpha
    tensor_nm2[1, 1] = a_nm * a_nm
    tensor_nm2[1, 2] = a_nm * a_nm * cos_alpha
    tensor_nm2[2, 0] = a_nm * a_nm * cos_alpha
    tensor_nm2[2, 1] = a_nm * a_nm * cos_alpha
    tensor_nm2[2, 2] = a_nm * a_nm

    return tensor_nm2


def gm_nm2(a_nm, b_nm, c_nm, beta_rad):
    tensor_nm2 = np.zeros((3, 3))
    cos_beta = np.cos(beta_rad)

    tensor_nm2[0, 0] = a_nm * a_nm
    tensor_nm2[1, 1] = b_nm * b_nm
    tensor_nm2[2, 2] = c_nm * c_nm
    tensor_nm2[0, 2] = a_nm * c_nm * cos_beta
    tensor_nm2[2, 0] = a_nm * c_nm * cos_beta

    return tensor_nm2


def ga_nm2(a_nm, b_nm, c_nm, alpha_rad, beta_rad, gamma_rad):
    tensor_nm2 = np.zeros((3, 3))

    tensor_nm2[0, 0] = a_nm * a_nm
    tensor_nm2[1, 1] = b_nm * b_nm
    tensor_nm2[2, 2] = c_nm * c_nm

    tensor_nm2[0, 1] = tensor_nm2[1, 0] = a_nm * b_nm * np.cos(gamma_rad)
    tensor_nm2[0, 2] = tensor_nm2[2, 0] = a_nm * c_nm * np.cos(beta_rad)
    tensor_nm2[1, 2] = tensor_nm2[2, 1] = b_nm * c_nm * np.cos(alpha_rad)

    return tensor_nm2
