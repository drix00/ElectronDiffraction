#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: crystal_system
   :synopsis: Crystallographic seven crystal systems.

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Crystallographic seven crystal systems.
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
from math import pi

# Third party modules.
import numpy as np

# Local modules.

# Project modules.

# Globals and constants variables.

class CrystalSystem():
    a_nm = None
    b_nm = None
    c_nm = None
    alpha_rad = None
    beta_rad = None
    gamma_rad = None

    def __init__(self, a_nm, b_nm, c_nm, alpha_rad, beta_rad, gamma_rad):
        self.a_nm = a_nm
        self.b_nm = b_nm
        self.c_nm = c_nm
        self.alpha_rad = alpha_rad
        self.beta_rad = beta_rad
        self.gamma_rad = gamma_rad

    def _compute_direct_metric_tensor(self):
        g_ij_nm2 = np.zeros((3, 3))

        g_ij_nm2[0, 0] = self.a_nm * self.a_nm
        g_ij_nm2[1, 1] = self.b_nm * self.b_nm
        g_ij_nm2[2, 2] = self.c_nm * self.c_nm

        g_ij_nm2[0, 1] = g_ij_nm2[1, 0] = self.a_nm * self.b_nm * np.cos(self.gamma_rad)
        g_ij_nm2[0, 2] = g_ij_nm2[2, 0] = self.a_nm * self.c_nm * np.cos(self.beta_rad)
        g_ij_nm2[1, 2] = g_ij_nm2[2, 1] = self.b_nm * self.c_nm * np.cos(self.alpha_rad)

        return g_ij_nm2

    @property
    def gij_nm2(self):
        return self._compute_direct_metric_tensor()


class Triclinic(CrystalSystem):
    pass


class Monoclinic(CrystalSystem):
    def __init__(self, a_nm, b_nm, c_nm, beta_rad):
        super().__init__(a_nm, b_nm, c_nm, pi/2.0, beta_rad, pi/2.0)


class Hexagonal(CrystalSystem):
    def __init__(self, a_nm, c_nm):
        super().__init__(a_nm, a_nm, c_nm, pi/2.0, pi/2.0, 2.0*pi/3.0)


class Rhombohedral(CrystalSystem):
    def __init__(self, a_nm, alpha_rad):
        super().__init__(a_nm, a_nm, a_nm, alpha_rad, alpha_rad, alpha_rad)


class Orthorhombic(CrystalSystem):
    def __init__(self, a_nm, b_nm, c_nm):
        super().__init__(a_nm, b_nm, c_nm, pi/2.0, pi/2.0, pi/2.0)


class Tetragonal(CrystalSystem):
    def __init__(self, a_nm, c_nm):
        super().__init__(a_nm, a_nm, c_nm, pi/2.0, pi/2.0, pi/2.0)


class Cubic(CrystalSystem):
    def __init__(self, a_nm):
        super().__init__(a_nm, a_nm, a_nm, pi/2.0, pi/2.0, pi/2.0)
