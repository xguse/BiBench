###################################################################################################
###---------------------------------------------------------------------------------------------###
### This file is part of bibench (Biclustering Benchmarking)			                        ###
### Copyright (c) 2011,                                                                         ###
### By:    Kemal Eren,                                                                          ###
##         Mehmet Deveci,                                                                       ###
###        Onur Kucuktunc,                                                                      ###
###        Umit V. Catalyurek                                                                   ###
###                                                                                             ###
###---------------------------------------------------------------------------------------------###
### For license info, please see the README.txt and LICENSE.txt files in the main directory.    ###
###---------------------------------------------------------------------------------------------###
###################################################################################################
"""
Command-line interface to BiBench's synthetic dataset creation
functionality.

"""

import sys
from bibench.algorithms import \
    bbc, biclust, cpb, opsm, fabia, qubic, coalesce, isa
from bibench.bicluster import write_biclusters
import argparse
from util import create_subparsers, get_kwargs
from numpy import genfromtxt


parser_parent = argparse.ArgumentParser(add_help=False,
                                     description='args common to all')
parser_parent.add_argument('datafile', action="store",
                           help='file containing data to bicluster')
parser_parent.add_argument('clusterfile', action="store",
                           help='file in which to store the found biclusters')

parser = argparse.ArgumentParser(
    description='Bicluster a dataset')
subparsers_adder = parser.add_subparsers(title='biclustering algorithms')

cmds = {
    'bbc': bbc.bbc,
    'cheng_church': biclust.cheng_church,
    'xmotifs': biclust.xmotifs,
    'bimax': biclust.bimax,
    'plaid': biclust.plaid,
    'spectral': biclust.spectral,
    'coalesce': coalesce.coalesce,
    'cpb': cpb.cpb,
    'fabia': fabia.fabia,
    'isa': isa.isa,
    'opsm': opsm.opsm,
    'qubic': qubic.qubic,
}

ignore = ['data']
create_subparsers(subparsers_adder.add_parser, cmds, [parser_parent], ignore=ignore)

if __name__ == '__main__':
    kwargs = get_kwargs(parser)
    datafile = kwargs.pop('datafile')
    data = genfromtxt(datafile)
    clusterfile = kwargs.pop('clusterfile')
    func = kwargs.pop('func')
    found = func(data, **kwargs)
    write_biclusters(found, clusterfile)
