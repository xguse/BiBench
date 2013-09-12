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
from bibench.datasets import synthetic
from bibench.bicluster import write_biclusters
from numpy import savetxt
import argparse
from util import create_subparsers, get_kwargs


parser_parent = argparse.ArgumentParser(add_help=False,
                                     description='args common to all')
parser_parent.add_argument('datafile', action="store",
                           help='file in which to store the generated dataset')
parser_parent.add_argument('clusterfile', action="store",
                           help='file in which to store the expected biclusters')
parser_parent.add_argument('-s', '--shuffle', action="store_true", default=False,
                           help='shuffle the dataset')

parser = argparse.ArgumentParser(
    description='Generate a dataset with implanted biclusters.')
subparsers_adder = parser.add_subparsers(title='data types')

cmds = {'const': synthetic.make_const_data,
        'shift': synthetic.make_shift_data,
        'scale': synthetic.make_scale_data,
        'shift_scale': synthetic.make_shift_scale_data,
        'plaid': synthetic.make_plaid_data,
        'fabia': synthetic.make_fabia_data,
        'isa': synthetic.make_isa_data}

# a list of arguments that can take multiple values
mult = ['bicluster_signals',
        'bicluster_noise']

ignore = ['dist',
          'shuffle']

create_subparsers(subparsers_adder.add_parser,
                  cmds,
                  [parser_parent],
                  ignore = ignore,
                  mult=mult)


if __name__ == '__main__':
    kwargs = get_kwargs(parser)
    datafile = kwargs.pop('datafile')
    clusterfile = kwargs.pop('clusterfile')
    func = kwargs.pop('func')
    data, expected = func(**kwargs)
    savetxt(datafile, data)
    write_biclusters(expected, clusterfile)
