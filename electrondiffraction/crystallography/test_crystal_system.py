#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. py:currentmodule:: test_crystal_system
   :synopsis: Tests for the module :py:mod:`crystal_system`

.. moduleauthor:: Hendrix Demers <hendrix.demers@mail.mcgill.ca>

Tests for the module :py:mod:`crystal_system`.
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
import electrondiffraction.crystallography.crystal_system as crystal_system

# Globals and constants variables.

class Test_crystal_system(unittest.TestCase):
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
        self.assert_(True)

    def test_gij_nm2(self):
        """
        Test the calculation of the direct metric tensor.

        From Marc De Graef Introduction to Conventional Transmission book (2003).
        """

        # Example 1.1
        g_ij_ref_nm2 = np.array([[0.25, 0.0, 0.0],
                                 [0.0, 0.25, 0.0],
                                 [0.0, 0.0, 1.0]])
        crystal = crystal_system.Tetragonal(0.5, 1.0)
        g_ij_nm2 = crystal.gij_nm2

        for g_ref_nm2, g_nm2 in zip(g_ij_ref_nm2.flat, g_ij_nm2.flat):
            self.assertAlmostEqual(g_ref_nm2, g_nm2, 5)

        #self.fail("Test if the testcase is working.")


if __name__ == '__main__':  # pragma: no cover
    import nose

    nose.runmodule()
