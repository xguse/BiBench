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
Biclustering algorithm wrappers for the R package 'biclust',
which includes several biclustering algorithms.

"""

import sys
import logging
from cStringIO import StringIO

from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm

from bibench.datasets.transform import is_discrete, is_binary

import rpy2.robjects as robjects
from rpy2.rinterface import RRuntimeError
import numpy

#enables automatic conversion from numpy to R:
import rpy2.robjects.numpy2ri
import pkg_resources
v = pkg_resources.get_distribution('rpy2').version
if v[0:3] >= '2.2':
    rpy2.robjects.numpy2ri.activate()


def _run_biclust_(function_name, data, **kwargs):
    """Convenience function for the various methods implemented in 'biclust'.

    Performs biclustering on the dataset and returns a set of biclusters.

    """
    #replace underscores with dots:
    keys = kwargs.keys()
    for key in keys:
        kwargs[key.replace("_", ".")] = kwargs.pop(key)

    robjects.r.library('biclust')

    #run biclustering
    biclust = robjects.r["biclust"]
    function = robjects.r[function_name]

    try:
        result = biclust(data, method=function_name, **kwargs)
    except RRuntimeError as e:
        logging.error(
            '{0} caught an R exception. Assuming no biclusters were found. Message: {1}'
            .format(function_name, e.message))
        return []

    #get rowXnumber array
    row_matrix = numpy.array(result.do_slot("RowxNumber"))

    #get numberXcolumn array
    col_matrix = numpy.array(result.do_slot("NumberxCol"))

    num_biclusters = row_matrix.shape[1]

    # a hack for Cheng and Church, which appears to sometimes get the transpose of
    # the column matrix
    if not num_biclusters == col_matrix.shape[0]:
        if num_biclusters == col_matrix.shape[1] and \
                row_matrix.shape[0] == data.shape[0] and \
                col_matrix.shape[0] == data.shape[1]:
            col_matrix = col_matrix.T
                
    if not num_biclusters == col_matrix.shape[0]:
        raise Exception(
            'There is a problem with the results returned by {0}'.format(function_name))

    #make list of biclusters
    biclusters = []
    for i in range(num_biclusters):
        rows_bools = row_matrix[:, i] != 0
        cols_bools = col_matrix[i, :] != 0

        rows = [index for index, elt in enumerate(rows_bools) if elt]
        cols = [index for index, elt in enumerate(cols_bools) if elt]

        biclusters.append(Bicluster(rows, cols, data=data))

    return biclusters


@bicluster_algorithm
def cheng_church(data, delta, alpha=1.5, number=100):
    """
    Cheng and Church algorithm.

    Args:
        * data: numpy.ndarray
        * delta: Maximum of accepted score.
        * alpha: Scaling factor.
        * number: Number of biclusters to find.

    Returns:
        A list of biclusters.
    """
    kwargs = locals()
    function_name = "BCCC"
    return _run_biclust_(function_name, **kwargs)


@bicluster_algorithm
def plaid(data,
          cluster="b",
          fit_model="y ~ m + a + b",
          background=True,
          row_release=0.7,
          col_release=0.7,
          shuffle=3,
          back_fit=0,
          max_layers=20,
          iter_startup=5,
          iter_layer=10,
          verbose=False):
    """
    The Plaid biclustering algorithm.

    Args:
        * data: numpy.ndarray
        * cluster: 'r', 'c' or 'b', to cluster rows, columns, or both.
        * fit_model: Formula to fit each layer.
            'm': const bicluster. 'a': const rows. 'b': const columns.
        * background: Consider a background layer present.
        * row_release: Threshold to prune rows in layers depending on homogeneity.
            Float in [0,1].
        * col_release: As row_release, but for columns. Float in [0,1].
        * shuffle: For computing statistical significance of a layer.
            Affects running time.
        * back_fit: Additional iterations for refining a layer.
        * max_layers: Maximum layers in the model.
        * iter_startup: Number of iterations to find starting values.
        * iter_layer: Number of iterations to find each layer.
        * verbose: if True, print progress.

    Returns:
        A list of biclusters.

    """
    kwargs = locals()
    kwargs["fit_model"] = robjects.r(fit_model)
    function_name = "BCPlaid"
    fnd = _run_biclust_(function_name, **kwargs)
    if len(fnd) == 1:
        if fnd[0].shape() == (1, 1):
            return []
    return fnd


@bicluster_algorithm
def bimax(data, minr=2, minc=2, number=100):
    """
    The BiMax biclustering algorithm. Searches for submatrices of ones
    in binary data.

    Notice: Bimax requres binary data. Method of binarization affects results.

    Args:
        * data: numpy.ndarray of zeroes and ones.
        * minr: Minimum rows in biclusters.
        * minc: Minimum columns in biclusters.
        * number: Number of biclusters to find.

    Returns:
        A list of biclusters.

    """
    kwargs = locals()
    if not is_binary(data):
        raise Exception('Bimax requires discrete data.')
    function_name = "BCBimax"
    return _run_biclust_(function_name, **kwargs)


@bicluster_algorithm
def spectral(data,
             normalization="log",
             numberOfEigenvalues=3,
             minr=2,
             minc=2,
             withinVar=1):
    """
    Finds checkerboard pattern in data using spectral decomposition

    Args:
        * data: numpy.ndarray
        * normalization: 'log', 'irrc', or 'bistochastization'.
        * numberOfEigenvalues: use the first eigenvalues.
            High numbers make runtime longer.
        * minr: minimum rows in a bicluster.
        * minc: minimum columns in a bicluster.
        * withinVar: biclusters with variance above this threshold are excluded.

    Returns:
        A list of biclusters.
    """
    kwargs = locals()
    function_name = "BCSpectral"

    stdout_save = sys.stdout
    sys.stdout = stdout = StringIO()
    found = _run_biclust_(function_name, **kwargs)
    sys.stdout = stdout_save
    stdout_val = stdout.getvalue()
    if stdout_val.find("No biclusters found") >= 0:
        return []
    return found


@bicluster_algorithm
def xmotifs(data, number=1, ns=200, nd=100, sd=5, alpha=0.05):
    """
    Finds biclusters with approximately constant rows/genes.

    Notice: XMotifs needs discrete data. The method of discretization
    affects the results.

    Args:
        * data: numpy.ndarray of ints.
        * number: number of biclusters to find
        * ns: number of seeds.
        * nd: number of determinants.
        * sd: size of discriminating set; generated for each seed.
        * alpha: scaling factor for column.

    Returns:
        A list of biclusters.

    """
    kwargs = locals()
    if not is_discrete(data):
        raise Exception('Xmotifs requires discrete data.')
    function_name = "BCXmotifs"
    return _run_biclust_(function_name, **kwargs)
