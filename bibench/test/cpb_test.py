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

class TestCpb(unittest.TestCase):

    #TODO: make this part of parent class
    def setUp(self):
        self.data, self.expected = bb.make_shift_scale_data()

    def test_cpb(self):
        nclus = 2
        result = bb.cpb(self.data, nclus, targetpcc=0.9)

        #cannot assume we will find exactly 'nclus' clusters, because
        #they are not returned by CPB if a row's PCC is below threshold.
        self.assertTrue(len(result) <= nclus)
        self.assertTrue(len(result) > 0)


    def test_cpb_filter(self):
        nclus = 100
        biclusters = bb.cpb(self.data, nclus, targetpcc=0.9)
        oldscore = bb.jaccard_list(self.expected, biclusters)

        filtered = bb.cpb_filter(biclusters, self.data, nclus, targetpcc=0.9)
        score = bb.jaccard_list(self.expected, filtered)

        self.assertTrue(len(filtered) <= len(biclusters))
        self.assertTrue(score >= oldscore)


if __name__ == "__main__":
    unittest.main()



