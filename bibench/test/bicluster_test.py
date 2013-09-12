####################################################################
###     ____  _ ____                  _                          ###
###    | __ )(_) __ )  ___ _ __   ___| |__                       ###
###    |  _ \| |  _ \ / _ \ '_ \ / __| '_ \                      ###
###    | |_) | | |_) |  __/ | | | (__| | | |                     ###
###    |____/|_|____/ \___|_| |_|\___|_| |_|                     ###
###                                                              ###
###--------------------------------------------------------------###
###                                                              ###
### This file is part of the BiBench package for biclustering    ###
### analysis.                                                    ###
###                                                              ###
### Copyright (c) 2011 by:                                       ###
###   * Kemal Eren,                                              ###
###   * Mehmet Deveci,                                           ###
###   * Umit V. Catalyurek                                       ###
###                                                              ###
###--------------------------------------------------------------###
###                                                              ###
### For license info, please see the README and LICENSE files    ###
### in the main directory.                                       ###
###                                                              ###
###--------------------------------------------------------------###

"""Unit tests for the 'bicluster' module"""

from __future__ import division

import unittest
import numpy as np

import bibench.all as bb
from bibench.bicluster import get_row_col_matrices, Bicluster

class BiclusterTest(unittest.TestCase):

    def test_get_bicluster(self):
        data = np.arange(60).reshape(10, 6)
        array = np.array([[25, 27, 28],
                          [37, 39, 40],
                          [55, 57, 58]])
        rows = (4, 6, 9)
        cols = (1, 3, 4)
        bicluster = Bicluster(rows, cols, data)
        self.assertTrue(np.alltrue(array == bicluster.array()))

    def test_bicluster_eq(self):
        bic_a = Bicluster([1, 2, 3], [1, 2, 3])
        bic_b = Bicluster([1, 2, 3], [1, 2, 3])
        self.assertEquals(bic_a, bic_b)

        data = np.arange(10)
        bic_b.data = data
        self.assertNotEquals(bic_a, bic_b)

        bic_a.data = data
        self.assertEquals(bic_a, bic_b)

        bic_b.data = np.arange(10)
        self.assertNotEquals(bic_a, bic_b)


    def test_get_row_col_matrices(self):
        exp_rows = np.vstack(np.array([1, 1, 0]))
        exp_cols = np.vstack(np.array([0, 1, 0]))

        data = np.random.randn(3, 3)
        biclusters = [Bicluster([0, 1], [1], data)]

        rowxnumber, colxnumber = get_row_col_matrices(biclusters)

        self.assertTrue((rowxnumber == exp_rows).all())
        self.assertTrue((colxnumber == exp_cols).all())


if __name__ == '__main__':
    unittest.main()
