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
Functions for data tranformations, such as binarization,
discretization, and normalization.

"""

import math
import numpy

import rpy2.robjects as robjects
import rpy2.interactive as r
from bibench.util import flatten

#enables automatic conversion from numpy to R:
import rpy2.robjects.numpy2ri
import pkg_resources
v = pkg_resources.get_distribution('rpy2').version
if v[0:3] >= '2.2':
    rpy2.robjects.numpy2ri.activate()


def _same_type_(data, orig):
    newdata = orig.copy()
    newdata[:,:] = data[:,:]
    return newdata


def _rfunction_(functionname, data, **kwargs):
    """
    get an R object for the data
    """
    r_data = robjects.Matrix(data)

    #get the function
    robjects.r.library('biclust')
    func = robjects.r[functionname]

    result = numpy.array(func(r_data, **kwargs))
    return _same_type_(result, data)


def is_discrete(data):
    """
    Checks if data is discrete.

    Args:
        * data: numpy.ndarray.

    """
    return data.dtype <= numpy.integer


def discretize(data, nof=10, quant=False, flip=True):
    """
    Discretizes the data matrix.

    Args:
        * data: numpy.ndarray.
        * nof: The number of levels in discretized matrix.
        * quant: Whether to use quantization
        * flip: If True, the largest level corresponds to the
            largest original value.

    Returns:
        numpy.ndarray; discretized version of 'data'.

    """
    kwargs = locals()
    kwargs.pop("flip")
    result = _rfunction_("discretize", **kwargs)
    if flip:
        result = result.max() - result
    return numpy.int8(result)


def is_binary(data):
    """
    Checks if data is binarized.

    Args:
        * data: numpy.ndarray

    Returns:
        bool

    """
    return data.dtype <= numpy.integer and \
        data.min() in (0, 1) and data.max() in (0, 1)


def binarize(data, threshold):
    """
    Binarizes the data matrix according to given threshold expression
    value.

    Args:
        * data: numpy.ndarray
        * threshold: The cutoff level for binarization.

    Returns:
        Binarized numpy.ndarray.

    """
    kwargs = locals()
    return numpy.int8(_rfunction_("binarize", **kwargs))


def _scoreatpercentile_(a, per):
    """
    My implementation, to stop using scipy.

    >>> _scoreatpercentile_(numpy.arange(100), 50)
    50.5

    """
    a = sorted(a)
    index = len(a) * per / 100.
    if index == int(math.floor(index)):
        index = int(math.floor(index))
        return (a[index] + a[index + 1]) / 2.
    return a[int(math.ceil(index))]


def binarize_quantile(data, quantile=0.5):
    """
    Binarizes the data matrix according to given quantile ratio.

    Args:
        * data: 2 dimensinal numpy.array format to represent data
            matrix that will be binarized.
        * quantile: the ratio of the 1's to total in the final
            binarized matrix. Default is 0.5

    Returns:
        Binarized numpy.ndarray.
    """
    thresh = _scoreatpercentile_(data.flatten(), quantile * 100)
    return binarize(data, thresh)


def densityOnes(data):
    """
    Returns the percentage of ones in the binary dataset.

    Args:
        * data: numpy.ndarray

    Returns:
        float

    """
    result = _rfunction_("densityOnes", data)
    return result[0]


def _standardize_helper_(data, axis):
    trans = lambda x: numpy.vstack(x) if axis == 1 else x
    return (data - trans(data.mean(axis=axis))) / trans(data.std(axis=axis))


def standardize(data):
    """
    Standardize entire dataset to have mean 0 and standard deviation 1.

    Args:
        * data: numpy.ndarray

    Returns:
        Standardized numpy.ndarray.

    """
    return _standardize_helper_(data, axis=None)


def standardize_cols(data):
    """
    Standardize columns to have mean 0 and standard deviation 1.

    Args:
        * data: numpy.ndarray.

    Returns:
        numpy.ndarray with standardized columns.

    """
    return _standardize_helper_(data, axis=0)


def standardize_rows(data):
    """
    Standardize rows to have mean 0 and standard deviation 1.

    Args:
        * data: numpy.ndarray.

    Returns:
        numpy.ndarray with standardized rows.

    """
    return _standardize_helper_(data, axis=1)


def log(data):
    """
    First shift so there are no zeros or negatives, then return the log
    of the data matrix, where log of each expression value is taken.

    Args:
        * data: numpy.ndarray
        * base: the base of the log function.

    Returns:
        numpy.ndarray.

    """
    if data.min() <= 0:
        data = data + abs(data.min()) + 1
    return numpy.log(data)


def _make_ranks_(length, nranks, down=False):
    """
    Returns a vector of ranks of length 'length', with
    each rank as equally represented as possible.

    >>> _make_ranks_(10, 3, down=False)
    [1, 1, 1, 1, 2, 2, 2, 3, 3, 3]

    >>> _make_ranks_(10, 3, down=True)
    [-1, -1, -1, -1, -2, -2, -2, -3, -3, -3]

    """
    if length < nranks:
        nranks = length #TODO: is this OK?
    assert nranks > 0
    if length < nranks:
        nranks = length
    neach, remain = divmod(length, nranks)
    ranks = [[i+1] * neach for i in range(nranks)]
    for rank in ranks[:remain]:
        rank.append(rank[0])
    ranks = flatten(ranks)
    if down:
        ranks = [-1 * r for r in ranks]
    return ranks


def _assign_ranks_(old, new, bools, nranks, down=False):
    """
    Modify the vector 'new' in-place, discretizing it to have
    'nranks' ranks.

    """
    regulated_vals = old[bools]
    if len(regulated_vals) == 0:
        return
    regulated_indices = numpy.arange(len(old))[bools]

    order = numpy.argsort(regulated_vals)[::-1]
    indices_rank_order = regulated_indices[order]

    ranks = _make_ranks_(len(order), nranks, down=down)
    for rank, index in zip(ranks, indices_rank_order):
        new[index] = rank


def qubic_discretize(data, quantile=0.06, nranks=1, up=True, down=True):
    """
    Quantize each row in the dataset according to qubic's method.

    Args:
        * data: numpy.ndarray.
        * quantile: Defines whether up or down regulated.
        * nranks: Number of discrete ranks in resulting dataset.
        * up: Whether to quantize upregulated elements.
        * down: Whether to quantize downregulated elements.

    Returns:
        numpy.ndarray.
    """
    new_data = numpy.zeros(data.shape, dtype=numpy.int32)
    ncols = data.shape[1]
    s = int(ncols * quantile) + 1
    for old, new in zip(data, new_data):
        median = numpy.median(old)
        lower = sorted(old)[s]
        upper = sorted(old)[ncols - s]
        d = min(median - lower, upper - median)

        if up:
            upregulated = old > median + d
            _assign_ranks_(old, new, upregulated, nranks, down=False)

        if down:
            downregulated = old < median - d
            _assign_ranks_(old, new, downregulated, nranks, down=True)
    return new_data


def qubic_binarize_up(data, quantile=0.06):
    """Binarize dataset so that upregulated elements are 1."""
    return qubic_discretize(data,
                            quantile=quantile,
                            nranks=1,
                            up=True,
                            down=False)


def qubic_binarize_down(data, quantile=0.06):
    """Binarize dataset so that downregulated elements are 1."""
    return -1 * qubic_discretize(data,
                                 quantile=quantile,
                                 nranks=1,
                                 up=False,
                                 down=True)


def remove_na_rows(data):
    """
    Returns the data with rows with all missing values removed.

    """
    notnan = numpy.isnan(data) != True
    keep = numpy.array([bool(numpy.any(r)) for r in notnan])
    return data[keep,:]


def pca_impute(data, method='bpca', scale='none', center=True, npcs=5):
    """
    Impute the missing data elements using PCA.

    Requires the 'pcaMethods' Bioconductor package.

    Args:
        * data
        * method: Which PCA algorithm to use. One of:
            * 'svd'
            * 'ppca'
            * 'bpca'
            * 'svdImpute'
            * 'nipals
            * 'robustPca'
        * scale:
        * center:
        * npcs: number of principal components.

    """
    kwargs = locals()
    r.importr('pcaMethods')
    data = remove_na_rows(data)
    r_data = robjects.Matrix(data)
    prepped = robjects.r['prep'](r_data, scale=scale, center=center)
    result = robjects.r['pca'](prepped, method=method, center=False, nPcs=npcs)
    imputed = numpy.array(robjects.r['completeObs'](result))
    return _same_type_(imputed, data)
