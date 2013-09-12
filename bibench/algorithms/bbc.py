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
Bayesian BiClustering (BBC) algorithm. Uses Gibbs sampling to find
biclusters fitting the Bayesian biclustering model.

"""
from bibench.algorithms.wrapper import \
    wrapper_helper, WrapperException
from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm
from bibench.datasets.io import write_expression_data
import numpy
import os
import subprocess

BINARY = 'BBC'

def bbc(data, nclus, norm_method='none', alpha=None):
    """
    Wrapper to the BBC binary.

    If 'nclus' is a list of integers, tries to determine the number of
    clusters in the dataset by performing multiple clusterings and
    choosing the one with the best BIC.

    sqrn sometimes causes a divide by zero if the data is too uniform.

    Args:
        * data: numpy.array; input data.

        * nclus: number of biclusters to find. Either an integer or
            a list.

        * norm_method: one of 'none', 'csn', 'rsn', 'irqn', 'sqrn'.

        * alpha: alpha% quartile used for IRQN or SQRN normalization.

    Returns: BiclusterList

    """
    """

    """
    kwargs = locals()
    try:
        results = [_bbc_(**dict(kwargs.items() + dict(nclus=k).items())) for k in nclus]
        idx = numpy.argmin([r.properties['bic'] for r in results])
        return results[idx]
    except TypeError:
        return _bbc_(**kwargs)


@bicluster_algorithm
def _bbc_(data, nclus, norm_method='none', alpha=None):
    kwargs = locals()
    assert norm_method in ['none', 'csn', 'rsn', 'iqrn', 'sqrn']
    if norm_method in ['iqrn', 'sqrn'] and alpha is None:
        raise WrapperException(
            "normalization method '{0}' requires alpha".format(norm_method))
    if not norm_method in ['iqrn', 'sqrn'] and alpha is not None:
        raise WrapperException(
            "alpha only used in quartile normalization: 'irqn' or 'sqrn'")

    result = wrapper_helper(BINARY,
                            _write_dataset_,
                            _read_results_,
                            _do_call_,
                            **kwargs)

    try:
        biclusters, props = result
        return BiclusterList(biclusters, properties=props)
    except ValueError:
        return BiclusterList([])


def _do_call_(data, datafile, results_dir, **kwargs):
    resultfile = os.path.join(results_dir, 'results.txt')
    command = "{binary} -i {0}" \
        " -k {nclus}" \
        " -o {1}" \
        " -n {norm_method}".format(datafile,
                                   resultfile,
                                   binary=BINARY,
                                   **kwargs)
    if kwargs['alpha'] is not None:
        command += " -r {alpha}".format(**kwargs)
    subprocess.check_call(command.split())


def _write_dataset_(data, filename):
    write_expression_data(data, filename, sep='\t')


def _read_result_file_(filename, data):
    biclusters = []
    with open(filename, 'r') as f:
        rows = []
        cols = []

        header = f.readline().split()
        properties = dict(nstable = int(header[5]),
                          likelihood = float(header[7]),
                          nparams = int(float(header[11])),
                          bic = float(header[13]))

        target = rows
        for line in f:
            if line[0:9] == "bicluster":
                if not line[9] == '1': #make sure we've read one biclustert
                    biclusters.append(Bicluster(rows, cols, data=data))
                    rows = []
                    cols = []
                f.next()
                continue
            elif line[0:3] == "row":
                target = rows
                continue
            elif line[0:3] == "col":
                target = cols
                continue
            else:
                v = int(line.split()[0]) - 1
                target.append(v)
        #ensure we get last bicluster.
        biclusters.append(Bicluster(rows, cols, data=data))
    return biclusters, properties


def _read_results_(results_dir, data):
    files = os.listdir(results_dir)
    assert len(files) is 1
    return _read_result_file_(os.path.join(results_dir, files[0]), data)
