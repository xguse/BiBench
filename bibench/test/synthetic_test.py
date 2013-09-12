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

import unittest
import numpy
import bibench.all as bb
import bibench.datasets.synthetic as synth

class SyntheticTest(unittest.TestCase):

    def test_isa(self):
        data, expected = bb.make_isa_data()
        self.assertTrue(len(expected) == 3)

    def test_fabia(self):
        noisy, expected = bb.make_fabia_data(1000, 100, 10, 5, 5, 5, 10, 3.0, 0.2, 2.0, 1.0, 0.2, 3.0, 1.0)
        self.assertTrue(len(expected) == 10)


    def test_shuffle(self):
        dataset = numpy.array([[1, 2, 3],
                               [4, 5, 6],
                               [7, 8, 9]])
        biclusters = [bb.Bicluster([0], [0, 1], dataset)]
        shuffled_data, shuffled_biclusters = synth._shuffle_(dataset, biclusters)
        b = shuffled_biclusters[0].array()[0]
        self.assertTrue(all(b == [1,2]) or all(b == [2,1]))


if __name__ == '__main__':
    unittest.main()

