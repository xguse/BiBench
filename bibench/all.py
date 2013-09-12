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

"""
A conveniance model to avoid tedious imports from various submodules.

Imports the most commonly used functions from all of the BiBench
project, making them available in one module.

"""

#import everything that is useful

from bibench.bicluster import \
    filter, Bicluster, BiclusterList, write_biclusters, read_biclusters

from bibench.util import bootstrap

from bibench.datasets.synthetic import \
    make_const_data, make_shift_data, make_scale_data, \
    make_shift_scale_data, make_fabia_data, make_isa_data, \
    make_plaid_data

from bibench.datasets.transform import \
    is_discrete, is_binary, discretize, binarize, \
    binarize_quantile, standardize, standardize_rows, \
    standardize_cols, log, qubic_discretize, qubic_binarize_up, \
    qubic_binarize_down, pca_impute

from bibench.datasets.io import read_expression_data, \
    write_expression_data, write_david_multilist, write_david_list, \
    write_bicoverlapper

from bibench.datasets.rdata import \
    get_r_data, get_bioc_data, get_gds_data, geo_query

from bibench.algorithms.bbc import bbc
from bibench.algorithms.biclust import \
    cheng_church, xmotifs, bimax, plaid, spectral
from bibench.algorithms.coalesce import coalesce
from bibench.algorithms.cpb import cpb, cpb_filter
from bibench.algorithms.fabia import \
    fabia, fabiap, fabias, mfsc, nmfdiv, nmfeu
from bibench.algorithms.isa import isa
from bibench.algorithms.opsm import opsm
from bibench.algorithms.qubic import qubic

from bibench.validation.external import \
    jaccard_list, prelic_list, f_measure_list, recovery_relevance_list

from bibench.validation.enrichment import enrichment, goid_annot

from bibench.visualization import \
    heatmap, parallel_coordinates, bubbleplot
