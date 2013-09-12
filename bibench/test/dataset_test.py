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

from __future__ import division

import unittest
import numpy as np

from bibench.bicluster import Bicluster

class DatasetTest(unittest.TestCase):

    def test_get_bicluster(self):
        data = np.arange(60).reshape(10, 6)
        array = np.array([[25, 27, 28],
                          [37, 39, 40],
                          [55, 57, 58]])
        rows = (4, 6, 9)
        cols = (1, 3, 4)
        bicluster = Bicluster(rows, cols, data)
        self.assertTrue(np.alltrue(array == bicluster.array()))

if __name__ == "__main__":
    unittest.main()
