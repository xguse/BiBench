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

class FabiaTest(unittest.TestCase):

    def test_fabia(self):
        nclus = 2
        data, expected = bb.make_fabia_data(100, 10, nclus, 5, 5, 5, 10, 3, .2, 2, 1, .2, 3, 1)
        results = bb.fabia(data, p=nclus)
        self.assertEquals(len(results), nclus)


if __name__ == '__main__':
    unittest.main()
