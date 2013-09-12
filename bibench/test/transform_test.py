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

import unittest

import bibench.all as bb
import numpy as np

class TestTransform(unittest.TestCase):
    def setUp(self):
        self.data = np.random.randn(20, 10)

    def test_binarize(self):
        result = bb.binarize(self.data, 0.2)
        self.assertTrue(bb.is_binary(result))

    def test_discretize(self):
        result = bb.discretize(self.data)
        self.assertTrue(bb.is_discrete(result))

    def test_qubic_discretize(self):
        result = bb.qubic_discretize(self.data)
        self.assertTrue(bb.is_discrete(result))

    def test_log(self):
        data = np.array([np.e ** 1, np.e ** 2, np.e ** 3])
        result = bb.log(data)
        self.assertTrue(np.all(result == np.array([1, 2, 3])))


if __name__ == "__main__":
    unittest.main()
