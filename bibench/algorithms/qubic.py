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
QUalitative BIClustering algorithm. Efficient algorithm for finding biclusters
with scaling patterns.
"""
from bibench.algorithms.wrapper import wrapper_helper
from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm
from bibench import util

import os.path
import re
import subprocess

number_regex = re.compile('[0-9]+')

BINARY = 'qubic'

@bicluster_algorithm
def qubic(data,
          nblocks=100,
          quantile=0.06,
          ranks=1,
          discrete=False,
          filtering=True,
          min_col_width=2,
          consistency_level=0.95):
    """
    QUBIC biclustering algorithm

    Args:
        * data: numpy.ndarray.
        * nblocks: Number of biclusters to report.
        * quantile: Quantile to use for discretization.
        * ranks: Number of ranks in discrete data.
        * discrete: True if the data is already discrete.
        * filtering: Whether to filter overlapping biclusters.
        * min_col_width: Minimum number of columns in a bicluster.
        * consistency_level: the minimum ratio between the number of
            identical valid symbols in a column and the total number
            of rows in the output

    Returns:
        A list of biclusters.

    """
    kwargs = locals()
    kwargs['filtering'] = 0 if filtering else 1
    return wrapper_helper(BINARY,
                          _write_dataset_,
                          _read_results_,
                          _do_call_,
                          **kwargs)


def _do_call_(data, datafile, results_dir, **kwargs):
    command = "{binary} -i {0}" \
        " -q {quantile}" \
        " -r {ranks}" \
        " -f {filtering}" \
        " -k {min_col_width}" \
        " -c {consistency_level}" \
        " -o {nblocks}".format(datafile, binary=BINARY, **kwargs)
    if kwargs['discrete']:
        command += ' -d'
    subprocess.check_call(command.split())

def _write_dataset_(data, filename):
    nrows, ncols = data.shape
    with open(filename, 'w') as f:
        #write first line
        line = ['o']
        line.extend(map(lambda x: 'cond{0}'.format(x), range(ncols)))
        f.write('\t'.join(line))
        f.write('\n')

        #write gene lines
        for i, line in enumerate(data):
            f.write('gene{0}\t'.format(i))
            f.write('\t'.join(map(str, line)))
            f.write('\n')

def _get_expected_(string, regex):
    matches = re.search(regex, string)
    number = number_regex.search(matches.group()).group()
    return int(number)

def _get_regex_(string_to_match):
    return re.compile('{0} \[[0-9]+\]:'.format(string_to_match),
                      flags=re.MULTILINE)

_gene_regex_ = _get_regex_("Genes")
_cond_regex_ = _get_regex_("Conds")

def _parse_bicluster_(string, gene_dict, cond_dict, data):
    expected_ngenes = _get_expected_(string, _gene_regex_)
    expected_nconds = _get_expected_(string, _cond_regex_)

    #split after the gene part
    after_genes = re.split(_gene_regex_, string)[1]

    #split into genes and conditions
    gene_lines, cond_lines = re.split(_cond_regex_, after_genes)

    cond_lines = cond_lines.split('\n')[0]

    rows = _handle_gene_lines_(gene_lines, gene_dict)
    cols = _handle_cond_lines_(cond_lines, cond_dict)

    assert len(rows) == expected_ngenes
    assert len(cols) == expected_nconds

    return Bicluster(rows, cols, data)


def _handle_gene_lines_( string, name_dict):
    return [name_dict[name] for name in string.split()]


def _handle_cond_lines_(string, name_dict):
    return [name_dict[name] for name in string.split()]


def _get_names_(results_dir):
    #map gene and condition names to rows/columns
    datafile = os.path.join(results_dir, 'data.txt')
    with open(datafile) as f:
        first = f.readline()
        rest = f.read()

        cond_names = first.split()[1:]
        gene_names = [line.split()[0] for line in rest.split('\n')
                      if line is not '']

        gene_dict = util.make_index_map(gene_names)
        cond_dict = util.make_index_map(cond_names)

        return gene_dict, cond_dict


def _read_results_(results_dir, data):
    results_dir = os.path.split(results_dir)[0]
    filename = os.path.join(results_dir, 'data.txt.blocks')
    with open(filename, 'r') as f:
        result_string = f.read()

    gene_dict, cond_dict = _get_names_(results_dir)

    start_regex = re.compile('^BC[0-9]+\s*S=[0-9]+$', flags=re.MULTILINE)
    bicluster_strings = re.split(start_regex, result_string)[1:]
    biclusters = [_parse_bicluster_(string, gene_dict, cond_dict, data)
                  for string in bicluster_strings]
    return biclusters
