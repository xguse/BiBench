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
FABIA biclustering algorithm.

"""

import numpy
from bibench import util
from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm
from rpy2 import robjects

@bicluster_algorithm
def fabia(data,
          p=5,
          alpha=0.1,
          cyc=500,
          spl=0.5,
          spz=0.5,
          random=1.0,
          center=2,
          norm=1,
          scale=0.0,
          lap=1.0):
    """
    Wrapper for the FABIA biclustering algorithm.

    Args:
        * data: numpy.ndarray.
        * p: number of hidden factors = number of biclusters.
        * alpha: sparseness loadings (0.1 - 1.0).
        * cyc: number of iterations.
        * spl: sparseness prior loadings (0.5 - 2.0) (Laplace).
        * spz: sparseness factors (0.5 - 2.0).
        * random: if <=0, then by SVD.
            if >0: random initialization of loadings in [-random, random].
        * center: data centering: 1 (mean), 2 (median), > 2 (mode), 0 (no).
        * norm: data normalization: 1 (0.75-0.25 quantile), >1 (var=1), 0 (no).
        * scale: loading vectors are scaled in each iteration to the given
            variance. 0.0 indicates non scaling.
        * lap: minimal value of the variational parameter.

    Returns:
        A list of biclusters.

    """
    params = locals()
    return _call_helper_('fabia', **params)


def fabiap(data,
           alpha=0.1,
           cyc=500,
           spl=0.5,
           spz=0.5,
           sL=0.6,
           sZ=0.6,
           random=1.0,
           center=2,
           norm=1,
           scale=0.0,
           lap=1.0):
    """Post-projection Fabia."""
    params = locals()
    return _call_helper_('fabiap', **params)


def fabias(data,
           p=5,
           alpha=0.6,
           cyc=500,
           spz=0.5,
           random=1.0,
           center=2,
           norm=1,
           lap=1.0):
    """Sparseness projection"""
    params = locals()
    return _call_helper_('fabias', **params)


def mfsc(data, p=5, cyc=500, sL=0.6, sZ=0.6, center=2, norm=1):
    """
    Sparse Matrix Factorization for bicluster analysis (MFSC)
    (Hochreiter et al., 2010).

    """
    params = locals()
    return _call_helper_('mfsc', **params)


def nmfdiv(data, p=5, cyc=100):
    """
    Non-negative Matrix Factorization with Kullaback-Leibler
    divergence as objective.

    """
    params = locals()
    return _call_helper_('nmfdiv', **params)


def nmfeu(data, p=5, cyc=100):
    """Non-negative Sparse Matrix Factorization with sparseness constraints."""
    params = locals()
    return _call_helper_('nmfeu', **params)


def _extract_biclusters_(fact, thresZ=0.5, thresL=None):
    params = dict()
    params['thresZ'] = thresZ
    if thresL is not None:
        params['thresL'] = thresL
    extract = robjects.r['extractBic']
    result = extract(fact, **params)

    data = result.rx('X')[0]
    numpy_data = numpy.array(data)
    row_dict = util.make_index_map(list(data.names[0]))
    col_dict = util.make_index_map(list(data.names[1]))

    # an R matrix; each row is a bicluster
    biclusters = []
    r_biclusters = result.rx('bic')[0]
    for b in range(1, r_biclusters.nrow + 1): #r matrices are 1-indexed
        entry = r_biclusters.rx(b, True)

        rownames = list(entry.rx('bixn')[0])
        colnames = list(entry.rx('biypn')[0])
        rows = [row_dict[r] for r in rownames]
        cols = [col_dict[c] for c in colnames]
        biclusters.append(Bicluster(rows, cols, numpy_data))
    return biclusters


def _call_helper_(function_name, data, **kwargs):
    params = kwargs
    params['X'] = data

    robjects.r.library('fabia')
    r_data = robjects.Matrix(data)
    func = robjects.r[function_name]
    factorization = func(**params)
    return _extract_biclusters_(factorization)
