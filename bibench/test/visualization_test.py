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
from bibench.bicluster import Bicluster
from bibench.visualization import _get_r_biclust_

class VisualizationTest(unittest.TestCase):

    #FIXME: always ends with Segmentation Fault or Illegal Instruction
    def test__get_r_biclust_(self):

        exp_rows = np.array([[1, 0, 1],
                             [1, 1, 0]], dtype=np.bool8)
        exp_cols = np.array([[1, 1, 0],
                             [1, 0, 1]], dtype=np.bool8)
        data = np.random.randn(2, 2)
        biclusters = [Bicluster([0, 1], [0, 1], data),
                      Bicluster([1], [0], data),
                      Bicluster([0], [1], data)]

        result = _get_r_biclust_(biclusters)

        rows = np.array(result.do_slot("RowxNumber"))
        cols = np.array(result.do_slot("NumberxCol"))
        cols = cols.T

        self.assertTrue((rows == exp_rows).all())
        self.assertTrue((cols == exp_cols).all())

if __name__ == '__main__':
    unittest.main()
