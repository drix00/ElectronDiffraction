#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: test_vector
   :synopsis: Tests for the module :py:mod:`vector`

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the module :py:mod:`vector`.
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
import unittest

# Third party modules.
import numpy as np

# Local modules.

# Project modules.
import electrondiffraction.crystallography.vector as vector
import electrondiffraction.crystallography.crystal_system as crystal_system

# Globals and constants variables.

class Test_vector(unittest.TestCase):
    """
    TestCase class for the module `${moduleName}`.
    """

    def setUp(self):
        """
        Setup method.
        """

        unittest.TestCase.setUp(self)

    def tearDown(self):
        """
        Teardown method.
        """

        unittest.TestCase.tearDown(self)

    def testSkeleton(self):
        """
        First test to check if the testcase is working with the testing framework.
        """

        #self.fail("Test if the testcase is working.")
        #self.assert_(True)

    def test_vector_dot_product(self):
        """
        Test the calculation of the direct metric tensor.

        From Marc De Graef Introduction to Conventional Transmission book (2003).
        """

        # Example 1.2
        vector_p = np.array([0.5, 0.0, 0.5])
        vector_q = np.array([0.5, 0.5, 0.0])
        crystal = crystal_system.Tetragonal(0.5, 1.0)
        magnitude_ref_nm = 5.0/16.0

        vector_d = vector_p - vector_q
        magnitude_nm = vector.dot_product(crystal, vector_d, vector_d)

        self.assertAlmostEqual(magnitude_ref_nm, magnitude_nm, 5)

        # Example 1.3
        vector_p = np.array([1.0, 2.0, 0.0])
        vector_q = np.array([3.0, 1.0, 1.0])
        crystal = crystal_system.Tetragonal(0.5, 1.0)
        magnitude_ref_nm = 5.0/4.0

        magnitude_nm = vector.dot_product(crystal, vector_p, vector_q)
        self.assertAlmostEqual(magnitude_ref_nm, magnitude_nm, 5)

        magnitude_nm = vector.dot_product(crystal, vector_q, vector_p)
        self.assertAlmostEqual(magnitude_ref_nm, magnitude_nm, 5)

        #self.fail("Test if the testcase is working.")

    def test_vector_distance(self):
        """
        Test the calculation of the direct metric tensor.

        From Marc De Graef Introduction to Conventional Transmission book (2003).
        """

        # Example 1.2
        vector_p = np.array([0.5, 0.0, 0.5])
        vector_q = np.array([0.5, 0.5, 0.0])
        crystal = crystal_system.Tetragonal(0.5, 1.0)
        d_ref_nm = np.sqrt(5.0)/4.0

        g_ij_nm2 = crystal.gij_nm2

        vector_d = vector_p - vector_q
        d_nm = vector.distance(crystal, vector_d, vector_d)

        self.assertAlmostEqual(d_ref_nm, d_nm, 5)

        #self.fail("Test if the testcase is working.")

    def test_points_distance(self):
        """
        Test the calculation of the direct metric tensor.

        From Marc De Graef Introduction to Conventional Transmission book (2003).
        """

        # Example 1.2
        point_p = np.array([0.5, 0.0, 0.5])
        point_q = np.array([0.5, 0.5, 0.0])
        crystal = crystal_system.Tetragonal(0.5, 1.0)
        d_ref_nm = np.sqrt(5.0)/4.0

        g_ij_nm2 = crystal.gij_nm2

        d_nm = vector.distance_points(crystal, point_p, point_q)

        self.assertAlmostEqual(d_ref_nm, d_nm, 5)

        #self.fail("Test if the testcase is working.")

    def test_vectors_angle(self):
        """
        Test the calculation of the direct metric tensor.

        From Marc De Graef Introduction to Conventional Transmission book (2003).
        """

        # Example 1.3
        vector_p = np.array([1.0, 2.0, 0.0])
        vector_q = np.array([3.0, 1.0, 1.0])
        crystal = crystal_system.Tetragonal(0.5, 1.0)
        angle_ref_deg = 53.300774799510123

        angle_rad = vector.angle_rad(crystal, vector_p, vector_q)
        angle_deg = np.degrees(angle_rad)
        self.assertAlmostEqual(angle_ref_deg, angle_deg, 6)

        angle_rad = vector.angle_rad(crystal, vector_q, vector_p)
        angle_deg = np.degrees(angle_rad)
        self.assertAlmostEqual(angle_ref_deg, angle_deg, 6)

        #self.fail("Test if the testcase is working.")

    def test_vectors_angle2(self):
        """
        Test the calculation of the direct metric tensor.

        From Marc De Graef Introduction to Conventional Transmission book (2003).
        """

        # Example 1.4
        vector_p = np.array([1.0, 2.0, 0.0])
        vector_q = np.array([3.0, 1.0, 1.0])
        crystal = crystal_system.Tetragonal(0.5, 1.0)
        angle_ref_deg = 53.300774799510123

        angle_rad = vector.angle2_rad(crystal, vector_p, vector_q)
        angle_deg = np.degrees(angle_rad)
        self.assertAlmostEqual(angle_ref_deg, angle_deg, 6)

        angle_rad = vector.angle2_rad(crystal, vector_q, vector_p)
        angle_deg = np.degrees(angle_rad)
        self.assertAlmostEqual(angle_ref_deg, angle_deg, 6)

        #self.fail("Test if the testcase is working.")


if __name__ == '__main__':  # pragma: no cover
    import nose

    nose.runmodule()
