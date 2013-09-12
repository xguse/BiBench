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

"""Unit tests for the 'validation' module"""

import unittest
import numpy
import bibench.all as bb
from bibench.bicluster import Bicluster

class TestValidation(unittest.TestCase):
    """
    Contains test cases for testing the validation functions in
    the 'validation' module.

    """
    data = numpy.random.randn(10, 10)
    list1 = [Bicluster([0, 1, 2, 3], [0, 1, 2, 3], data)]
    list2 = [Bicluster([2, 3, 4, 5], [2, 3, 4, 5], data)]

    def test_prelic(self):
        rel, rec = bb.prelic_list(self.list1, self.list1)
        self.assertEqual(rel, 1)
        self.assertEqual(rec, 1)

        rel, rec = bb.prelic_list(self.list1, self.list2)
        self.assertAlmostEqual(rel, 1/3)
        self.assertAlmostEqual(rec, 1/3)


    def test_fmeasure(self):
        rel, rec = bb.f_measure_list(self.list1, self.list1, modified=False)
        self.assertEquals(rel, 1)
        self.assertEquals(rec, 1)

        sens = 4 / 16
        spec = (100 - 28) / (100 - 16)

        expected = 2 * (sens * spec) / (sens + spec)

        rel, rec = bb.f_measure_list(self.list1, self.list2, modified=False)
        self.assertAlmostEqual(rel, expected)
        self.assertAlmostEqual(rec, expected)



    def test_modified_fmeasure(self):
        rel, rec = bb.f_measure_list(self.list1, self.list1, modified=True)
        self.assertEquals(rel, 1)
        self.assertEquals(rec, 1)


    def test_bicluster_jaccard(self):
        rel, rec = bb.jaccard_list(self.list1, self.list1)
        self.assertEquals(rel, 1)
        self.assertEquals(rec, 1)

        expected = 4 / 28

        rel, rec = bb.jaccard_list(self.list1, self.list2)
        self.assertAlmostEqual(rel, expected)
        self.assertAlmostEqual(rec, expected)


    def test_recovery_and_relevance(self):
        rel, rec = bb.recovery_relevance_list(self.list1, self.list1)
        self.assertEquals(rel, 1)
        self.assertEquals(rec, 1)

        rel, rec = bb.recovery_relevance_list(self.list1, self.list2)
        self.assertAlmostEqual(rel, 0.25)
        self.assertAlmostEqual(rec, 0.25)


if __name__ == "__main__":
    unittest.main()
