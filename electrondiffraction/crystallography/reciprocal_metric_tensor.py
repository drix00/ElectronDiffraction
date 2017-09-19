#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: reciprocal_metric_tensor
   :synopsis: Equations to compute the reciprocal metric tensor for each crystal system.

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Equations to compute the reciprocal metric tensor for each crystal system.
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


def grc_1_nm2(a_nm):
    tensor_1_nm2 = np.zeros((3, 3))

    tensor_1_nm2[0, 0] = 1.0 / (a_nm * a_nm)
    tensor_1_nm2[1, 1] = 1.0 / (a_nm * a_nm)
    tensor_1_nm2[2, 2] = 1.0 / (a_nm * a_nm)

    return tensor_1_nm2


def grt_1_nm2(a_nm, c_nm):
    tensor_1_nm2 = np.zeros((3, 3))

    tensor_1_nm2[0, 0] = 1.0 / (a_nm * a_nm)
    tensor_1_nm2[1, 1] = 1.0 / (a_nm * a_nm)
    tensor_1_nm2[2, 2] = 1.0 / (c_nm * c_nm)

    return tensor_1_nm2


def gro_1_nm2(a_nm, b_nm, c_nm):
    tensor_1_nm2 = np.zeros((3, 3))

    tensor_1_nm2[0, 0] = 1.0 / (a_nm * a_nm)
    tensor_1_nm2[1, 1] = 1.0 / (b_nm * b_nm)
    tensor_1_nm2[2, 2] = 1.0 / (c_nm * c_nm)

    return tensor_1_nm2


def grh_1_nm2(a_nm, c_nm):
    tensor_1_nm2 = np.zeros((3, 3))

    tensor_1_nm2[0, 0] = 4.0 / (3.0 * a_nm * a_nm)
    tensor_1_nm2[1, 1] = 4.0 / (3.0 * a_nm * a_nm)
    tensor_1_nm2[2, 2] = 1.0 / (c_nm * c_nm)
    tensor_1_nm2[0, 1] = 2.0 / (3.0 * a_nm * a_nm)
    tensor_1_nm2[1, 0] = 2.0 / (3.0 * a_nm * a_nm)

    return tensor_1_nm2


def grr_1_nm2(a_nm, alpha_rad):
    tensor_1_nm2 = np.zeros((3, 3))
    cos_alpha = np.cos(alpha_rad)
    factor_tan2_alpha_2 = (1.0 - np.tan(alpha_rad / 2.0) * np.tan(alpha_rad / 2.0)) / 2.0

    tensor_1_nm2[0, 0] = 1.0 + cos_alpha
    tensor_1_nm2[0, 1] = -factor_tan2_alpha_2
    tensor_1_nm2[0, 2] = -factor_tan2_alpha_2
    tensor_1_nm2[1, 0] = -factor_tan2_alpha_2
    tensor_1_nm2[1, 1] = 1.0 + cos_alpha
    tensor_1_nm2[1, 2] = -factor_tan2_alpha_2
    tensor_1_nm2[2, 0] = -factor_tan2_alpha_2
    tensor_1_nm2[2, 1] = -factor_tan2_alpha_2
    tensor_1_nm2[2, 2] = 1.0 + cos_alpha

    W = _compute_W(a_nm, cos_alpha)

    tensor_1_nm2 = tensor_1_nm2 / (W * W)

    return tensor_1_nm2


def _compute_W(a_nm, cos_alpha):
    W = a_nm * a_nm * (1.0 + cos_alpha - 2.0 * cos_alpha * cos_alpha)
    return W


def grm_1_nm2(a_nm, b_nm, c_nm, beta_rad):
    tensor_1_nm2 = np.zeros((3, 3))
    cos_beta = np.cos(beta_rad)
    sin2_beta = np.sin(beta_rad) * np.sin(beta_rad)

    tensor_1_nm2[0, 0] = 1.0 / (a_nm * a_nm * sin2_beta)
    tensor_1_nm2[1, 1] = 1.0 / (b_nm * b_nm)
    tensor_1_nm2[2, 2] = 1.0 / (c_nm * c_nm * sin2_beta)
    tensor_1_nm2[0, 2] = -cos_beta / (a_nm * c_nm * sin2_beta)
    tensor_1_nm2[2, 0] = -cos_beta / (a_nm * c_nm * sin2_beta)

    return tensor_1_nm2


def gra_1_nm2(a_nm, b_nm, c_nm, alpha_rad, beta_rad, gamma_rad):
    tensor_1_nm2 = np.zeros((3, 3))

    tensor_1_nm2[0, 0] = b_nm * b_nm * c_nm * c_nm * np.sin(alpha_rad) * np.sin(alpha_rad)
    tensor_1_nm2[1, 1] = a_nm * a_nm * c_nm * c_nm * np.sin(beta_rad) * np.sin(beta_rad)
    tensor_1_nm2[2, 2] = a_nm * a_nm * b_nm * b_nm * np.sin(gamma_rad) * np.sin(gamma_rad)

    tensor_1_nm2[0, 1] = tensor_1_nm2[1, 0] = a_nm * b_nm * c_nm * c_nm * _compute_F(alpha_rad, beta_rad, gamma_rad)
    tensor_1_nm2[0, 2] = tensor_1_nm2[2, 0] = a_nm * b_nm * b_nm * c_nm * _compute_F(gamma_rad, alpha_rad, beta_rad)
    tensor_1_nm2[1, 2] = tensor_1_nm2[2, 1] = a_nm * a_nm * b_nm * c_nm * _compute_F(beta_rad, gamma_rad, alpha_rad)

    omega = _compute_omega(a_nm, b_nm, c_nm, alpha_rad, beta_rad, gamma_rad)

    tensor_1_nm2 = tensor_1_nm2 / (omega * omega)
    return tensor_1_nm2


def _compute_F(a_rad, b_rad, g_rad):
    Fabg = np.cos(a_rad) * np.cos(b_rad) - np.cos(g_rad)
    return Fabg


def _compute_omega(a_nm, b_nm, c_nm, alpha_rad, beta_rad, gamma_rad):
    factor = a_nm * b_nm * c_nm
    factor = factor * factor
    term1 = np.cos(alpha_rad) * np.cos(alpha_rad) + np.cos(beta_rad) * np.cos(beta_rad) + \
            np.cos(gamma_rad) * np.cos(gamma_rad)
    term2 = 2.0 * np.cos(alpha_rad) * np.cos(beta_rad) * np.cos(gamma_rad)

    omega = factor * (1.0 - term1 + term2)
    return omega
