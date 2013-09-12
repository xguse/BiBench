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

class TestBiclust(unittest.TestCase):

    def setUp(self):
        self.data = bb.get_r_data('SyntrenEcoli', 'biclust')

    def test_cc(self):
        result = bb.cheng_church(self.data, delta=3, number=1)
        self.assertEquals(len(result), 1)

    def test_cc_background(self):
        """
        Added after finding a ?bug? in which Cheng and Church's
        column matrix would be transposed for certain datasets, and
        not for others.

        It seems to happen when Cheng and Church finds the whole
        dataset as a bicluster; this often happens when delta is large.

        """
        data, exp = bb.make_shift_data(background_scale=0)
        fnd = bb.cheng_church(self.data, delta=0.5)
        

    def test_plaid(self):
        result = bb.plaid(self.data)
        self.assertTrue(len(result) > 0)


    def test_bimax(self):
        result = bb.bimax(bb.binarize(self.data, 1))
        self.assertTrue(len(result) > 0)

    def test_spectral(self):
        result = bb.spectral(self.data)
        self.assertTrue(len(result) > 0)

    #FIXME: this test fails nondeterministically
    def test_xmotifs(self):
        result = bb.xmotifs(bb.discretize(self.data))
        self.assertTrue(len(result) > 0)


if __name__ == "__main__":
    unittest.main()
