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
Methods for creating synthetic datasets with implanted biclusters.

"""

from bibench.bicluster import Bicluster

import numpy

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
import pkg_resources
v = pkg_resources.get_distribution('rpy2').version
if v[0:3] >= '2.2':
    rpy2.robjects.numpy2ri.activate()

import random


def improper_normal(loc=0, scale=1, size=None):
    """
    Same as numpy.random.normal, but if scale is 0, returns just the
    mean 'loc'.

    """
    result = None
    if scale == 0:
        if size is None:
            result = loc
        else:
            result = numpy.zeros(size) + loc
    else:
        result = numpy.random.normal(loc, scale, size)
    return result


def _add_bicluster_noise_(data, rowmatrix, colmatrix, stdevs):
    """Adds noise from normal(0, stdevs[i]) to the bicluster i in the data."""
    noisy = numpy.copy(data)
    for rows, cols, scale in zip(rowmatrix.T, colmatrix.T, stdevs):
        if scale > 0:
            shape = noisy[rows > 0][:, cols > 0].shape
            noisy[numpy.outer(rows, cols) > 0] += \
                numpy.random.normal(scale=scale, size=shape).flatten()
    return noisy


def add_noise(data, scale):
    """Adds normal noise with 0 mean and 'scale' standard deviation."""
    return data + improper_normal(scale=scale, size=data.shape)


def _shuffle_(data, expected, new_rows=None, new_cols=None):
    """
    Shuffles the dataset while preserving biclusters.

    Args:
        * data: numpy.ndarray
        * expected: list of biclusters.
        * new_rows: Shuffled row indices; if None, randomly generated.
        * new_cols: Shuffled column indices; if None, randomly generated.

    Returns:
        The tuple (shuffled_data, shuffled_biclusters) where shuffled_data
        is a shuffled version of the input dataset, and shuffled_biclusters
        is a list of biclusters corresponding to the new biclusters in
        the shuffled dataset.

    """
    nrows, ncols = data.shape
    if new_rows is None:
        new_rows = range(nrows)
        random.shuffle(new_rows)
    if new_cols is None:
        new_cols = range(ncols)
        random.shuffle(new_cols)

    shuffled_data = data[new_rows].T[new_cols].T
    shuffled_biclusters = []
    for b in expected:
        new_b_rows = [new_rows.index(r) for r in b.rows]
        new_b_cols = [new_cols.index(c) for c in b.cols]
        shuffled_biclusters.append(Bicluster(new_b_rows,
                                             new_b_cols,
                                             shuffled_data))
    return shuffled_data, shuffled_biclusters



def _make_row_matrix_(nrows, nclusts, nclust_rows, noverlap_rows):
    """
    Make a matrix with rows representing the rows of the dataset, and
    columns representing biclusters.  If matrix[i, j] == 1, then row i
    is in bicluster j.

    """
    #check sanity of arguments
    if nclust_rows + (nclust_rows - noverlap_rows) * (nclusts - 1) > nrows:
        raise Exception('biclusters are too large to fit in the dataset.')

    matrix = numpy.zeros((nrows, nclusts))
    for i in range(nclusts):
        start = (nclust_rows * i) - (noverlap_rows * i)
        stop = start + nclust_rows
        matrix[start:stop, i] = 1
    return matrix


def _make_col_matrix_(ncols, nclusts, nclust_cols, noverlap_cols):
    """
    Make a matrix with rows representing columns, and columns
    representing biclusters.  If matrix[i, j] == 1, then column i is
    in bicluster j.

    """
    return _make_row_matrix_(ncols, nclusts, nclust_cols, noverlap_cols)


def _make_expected_biclusters_(row_matrix, col_matrix, data):
    """
    Given the output of _make_row_matrix_() and _make_col_matrix_(),
    make a list of Biclusters.

    """
    nclust = row_matrix.shape[1]
    assert nclust == col_matrix.shape[1]

    biclusters = []
    for row_line, col_line in zip(row_matrix.T, col_matrix.T):
        rows = list(numpy.where(row_line > 0)[0])
        cols = list(numpy.where(col_line > 0)[0])
        biclusters.append(Bicluster(rows, cols, data))
    return biclusters


def _make_biclusters_(row_matrix,
                      col_matrix,
                      bicluster_colbase,
                      bicluster_rowshift,
                      bicluster_rowscale):
    nrows, nclusts = row_matrix.shape
    ncols, nclusts_check = col_matrix.shape

    data = numpy.zeros((nrows, ncols))
    for row, col, base, shift, scale in zip(row_matrix.T,
                                            col_matrix.T,
                                            bicluster_colbase.T,
                                            bicluster_rowshift.T,
                                            bicluster_rowscale.T):
        bicluster = numpy.outer(scale, base) + numpy.vstack(shift)
        bool_matrix = numpy.outer(row, col) > 0

        #hack for True matrix
        data[bool_matrix] = bicluster[bicluster == bicluster]
    return data


def _do_both_overlap_(noverlap_rows, noverlap_cols):
    return noverlap_rows > 0 and noverlap_cols > 0


def _general_model_(nrows, ncols,
                    nclusts, nclustrows, nclustcols,
                    noverlap_rows, noverlap_cols,
                    colbase, rowshift, rowscale,
                    background_generator):
    """
    Make a dataset, given appropriate parameters.

    If rows or columns do not overlap: colbase should be ncolsx1 and
    rowshift, rowscale should be nrowsx1.  Else if both overlap:
    colbase should be nclustcolsxnclusts and rowshift, rowscale should
    be nclustrowsx1.

    """

    k = numpy.identity(nclusts)
    row_matrix = _make_row_matrix_(nrows, nclusts, nclustrows, noverlap_rows)
    col_matrix = _make_col_matrix_(ncols, nclusts, nclustcols, noverlap_cols)

    data = numpy.dot(numpy.dot(row_matrix, k), col_matrix.T)
    bicluster_bools = data > 0
    background_bools = data == 0

    background = background_generator(nrows, ncols)

    if _do_both_overlap_(noverlap_rows, noverlap_cols):
         #FIXME: check exactly what the dimensions are.
        biclusters = numpy.outer(rowscale, colbase) + numpy.vstack(rowshift)
    else:
        biclusters = _make_biclusters_(row_matrix,
                                       col_matrix,
                                       colbase,
                                       rowshift,
                                       rowscale)

    data[bicluster_bools] = biclusters[bicluster_bools]
    data[background_bools] = background[background_bools]

    return data, row_matrix, col_matrix


def _set_defaults_(nrows,
                   ncols,
                   nclusts,
                   nclustrows,
                   nclustcols,
                   bicluster_noise):
    if nclustrows is None:
        nclustrows = divmod(nrows, nclusts)[0]
        if nclustrows == nrows:
            nclustrows = int(nrows / 2)
    if nclustcols is None:
        nclustcols = divmod(ncols, nclusts)[0]
        if nclustcols == ncols:
            nclustcols = int(ncols / 2)
    if bicluster_noise is None:
        bicluster_noise = [0] * nclusts
    return nclustrows, nclustcols, bicluster_noise


def _make_data_(nrows, ncols, nclusts,
                noverlap_rows, noverlap_cols,
                colbase, rowshift, rowscale,
                noise, nclustrows=None, nclustcols=None, bicluster_noise=None,
                background_loc=0, background_scale=1, shuffle=False,
                dist=improper_normal):
    """Takes the matrices, applies the model, and adds noise."""
    background_generator = lambda x, y: dist(loc=background_loc,
                                             scale=background_scale,
                                             size=(x, y))

    data, row_matrix, col_matrix = _general_model_(nrows,
                                                   ncols,
                                                   nclusts,
                                                   nclustrows,
                                                   nclustcols,
                                                   noverlap_rows,
                                                   noverlap_cols,
                                                   colbase,
                                                   rowshift,
                                                   rowscale,
                                                   background_generator)

    if noise > 0:
        data = add_noise(data, scale=noise)

    #add extra noise to biclusters
    data = _add_bicluster_noise_(data, row_matrix, col_matrix, bicluster_noise)

    expected = _make_expected_biclusters_(row_matrix, col_matrix, data)

    if shuffle:
        data, expected = _shuffle_(data, expected)

    return data, expected



def make_const_data(nrows=300, ncols=50,
                    nclusts=3, nclustrows=None, nclustcols=None,
                    noverlap_rows=0, noverlap_cols=0,
                    background_loc=0, background_scale=1,
                    bicluster_signals=None,
                    noise=0, bicluster_noise=None,
                    shuffle=False,
                    dist=improper_normal):
    """
    Create a synthetic dataset with constant biclusters.

    Args: see make_shift_scale_data()
        * bicluster_signals: Constant bicluster expression values. May
          be either an int or a list of ints, one for each
          bicluster. Defaults to [0] * nclusts

    Returns:
        (numpy.ndarray, list of expected biclusters)

    """
    nclustrows, nclustcols, bicluster_noise = _set_defaults_(nrows,
                                                             ncols,
                                                             nclusts,
                                                             nclustrows,
                                                             nclustcols,
                                                             bicluster_noise)
    if bicluster_signals is None:
        bicluster_signals = [0] * nclusts
    if type(bicluster_signals) is int:
        bicluster_signals = [bicluster_signals] * nclusts

    assert len(bicluster_signals) == nclusts

    if _do_both_overlap_(noverlap_rows, noverlap_cols):
        signal = bicluster_signals[0]
        assert all([i == signal for i in bicluster_signals])
        colbase = numpy.array([0] * ncols)
        indexlim = nclustcols + (nclustcols - noverlap_cols) * (nclusts - 1)
        colbase[0:indexlim] = signal
        for i, signal in enumerate(bicluster_signals):
            colbase[i * nclustcols : (i + 1) * nclustcols] = signal
        rowshift = numpy.array([0] * nrows)
        rowscale = numpy.array([1] * nrows)
    else: #use seperate bases, shifting, and scaling parameters
        colbase = numpy.zeros((nclustcols, nclusts))
        for i, signal in enumerate(bicluster_signals):
            colbase[:,i] += signal
        rowshift = numpy.zeros((nclustrows, nclusts))
        rowscale = numpy.ones((nclustrows, nclusts))

    return _make_data_(nrows, ncols, nclusts,
                       noverlap_rows, noverlap_cols,
                       colbase, rowshift, rowscale,
                       noise, nclustrows, nclustcols, bicluster_noise,
                       background_loc, background_scale, shuffle, dist)


def make_shift_data(nrows=300, ncols=50,
                    nclusts=3, nclustrows=None, nclustcols=None,
                    noverlap_rows=0, noverlap_cols=0,
                    background_loc=0, background_scale=1,
                    base_loc=0, base_scale=1,
                    shift_loc=0, shift_scale=1,
                    noise=0, bicluster_noise=None,
                    shuffle=False,
                    dist=improper_normal):


    """
    Create a synthetic dataset with additive(shift) biclusters.

    Args: see make_shift_scale_data()

    """
    nclustrows, nclustcols, bicluster_noise = _set_defaults_(nrows,
                                                             ncols,
                                                             nclusts,
                                                             nclustrows,
                                                             nclustcols,
                                                             bicluster_noise)
    if _do_both_overlap_(noverlap_rows, noverlap_cols):
        colbase = dist(loc=base_loc, scale=base_scale, size=ncols)
        rowshift = dist(loc=shift_loc, scale=shift_scale, size=nrows)
        rowscale = numpy.array([1] * nrows)
    else:
        colbase = dist(loc=base_loc,
                       scale=base_scale,
                       size=(nclustcols, nclusts))
        rowshift = dist(loc=shift_loc,
                        scale=shift_scale,
                        size=(nclustrows, nclusts))
        rowscale = numpy.ones((nclustrows, nclusts))

    return _make_data_(nrows, ncols, nclusts,
                       noverlap_rows, noverlap_cols,
                       colbase, rowshift, rowscale,
                       noise, nclustrows, nclustcols, bicluster_noise,
                       background_loc, background_scale, shuffle,
                       dist)


def make_scale_data(nrows=300, ncols=50,
                    nclusts=3, nclustrows=None, nclustcols=None,
                    noverlap_rows=0, noverlap_cols=0,
                    background_loc=0, background_scale=1,
                    base_loc=0, base_scale=1,
                    scale_loc=0, scale_scale=1,
                    noise=0, bicluster_noise=None, shuffle=False,
                    dist=improper_normal):
    """
    Create a synthetic dataset with multiplicative (scale)
    biclusters.

    Args: see make_shift_scale_data()

    """
    nclustrows, nclustcols, bicluster_noise = _set_defaults_(nrows,
                                                             ncols,
                                                             nclusts,
                                                             nclustrows,
                                                             nclustcols,
                                                             bicluster_noise)

    if _do_both_overlap_(noverlap_rows, noverlap_cols):
        colbase = dist(loc=base_loc, scale=base_scale, size=ncols)
        rowshift = numpy.array([0] * nrows)
        rowscale = dist(loc=scale_loc, scale=scale_scale, size=nrows)
    else:
        colbase = dist(loc=base_loc,
                       scale=base_scale,
                       size=(nclustcols,
                             nclusts))
        rowshift = numpy.zeros((nclustrows, nclusts))
        rowscale = dist(loc=scale_loc,
                        scale=scale_scale,
                        size=(nclustrows,
                              nclusts))

    return _make_data_(nrows, ncols, nclusts,
                       noverlap_rows, noverlap_cols,
                       colbase, rowshift, rowscale,
                       noise, nclustrows, nclustcols, bicluster_noise,
                       background_loc, background_scale, shuffle,
                       dist)


def make_shift_scale_data(nrows=300, ncols=50,
                          nclusts=3, nclustrows=None, nclustcols=None,
                          noverlap_rows=0, noverlap_cols=0,
                          background_loc=0, background_scale=1,
                          base_loc=0, base_scale=1,
                          shift_loc=0, shift_scale=1,
                          scale_loc=0, scale_scale=1,
                          noise=0, bicluster_noise=None,
                          shuffle=False,
                          dist=improper_normal):
    """
    Create a synthetic dataset with shift-scale biclusters.

    Args:
        * nrows: Number of rows in dataset.
        * ncols: Number of columns in dataset.
        * nclusts: Number of biclusters.
        * nclustrows: Number of bicluster rows.
            Defaults to nrows / nclusts, rounded down
        * nclustcols: Number of bicluster columns.
            Defaults to ncols / nclusts, rounded down
        * noverlap_rows: Number of overlapping rows for each bicluster.
        * noverlap_cols: Number of overlapping columns for each bicluster.
        * background_loc, background_scale: Location and scale parameters for
            drawing background values from 'dist'.
        * base_loc, base_scale: Parameters for drawing base values.
        * shift_loc, shift_scale: Parameters for drawing shifting values.
        * scale_loc, scale_scale: Parameters for drawing scaling values.
        * noise: Scale parameter for distribution from which noise is drawn.
        * bicluster_noise: same as noise, but applied only to each bicluster.
            Defaults to [0] * nclusts
        * shuffle: If true, shuffle the rows and columns of resulting dataset.
        * dist: A distribution function that takes arguments (loc, scale, size).
            Used for generating random vectors and matrices.

    """
    nclustrows, nclustcols, bicluster_noise = _set_defaults_(nrows,
                                                             ncols,
                                                             nclusts,
                                                             nclustrows,
                                                             nclustcols,
                                                             bicluster_noise)
    if _do_both_overlap_(noverlap_rows, noverlap_cols):
        colbase = dist(loc=base_loc, scale=base_scale, size=ncols)
        rowshift = dist(loc=shift_loc, scale=shift_scale, size=nrows)
        rowscale = dist(loc=scale_loc, scale=scale_scale, size=nrows)
    else:
        colbase = dist(loc=base_loc,
                       scale=base_scale,
                       size=(nclustcols, nclusts))
        rowshift = dist(loc=shift_loc,
                        scale=shift_scale,
                        size=(nclustrows, nclusts))
        rowscale = dist(loc=scale_loc,
                        scale=scale_scale,
                        size=(nclustrows, nclusts))

    return _make_data_(nrows, ncols, nclusts,
                       noverlap_rows, noverlap_cols,
                       colbase, rowshift, rowscale,
                       noise, nclustrows, nclustcols, bicluster_noise,
                       background_loc, background_scale, shuffle,
                       dist)


def make_plaid_data(nrows=300, ncols=50,
                    nclusts=3, nclustrows=None, nclustcols=None,
                    noverlap_rows=0, noverlap_cols=0,
                    default_scale = 0.5, #used when scale not given
                    background_loc = 0,
                    background_scale = 0.5,
                    cluster_locs = None,
                    cluster_scales = None,
                    row_locs = None,
                    row_scales = None,
                    col_locs = None,
                    col_scales = None,
                    error_scales = None, #defaults to 0
                    noise = 0,
                    shuffle=False,
                    dist=improper_normal):
    """
    Create a synthetic dataset with plaid biclusters.

    Args: see make_shift_scale_data()
        * default_scale: plaid default scale; used for any unspecified scale
            parameter.
        * cluster_locs, cluster_scales: Location and scale parameters used to
            generate each cluster effect.
        * row_locs, row_scales: Parameters used to generate row effects.
        * col_locs, col_scales: Parameters used to generate column effects.
        * error_scales: Scaling parameters used for cluster error.

    """

    nclustrows, nclustcols, bicluster_noise = _set_defaults_(nrows,
                                                             ncols,
                                                             nclusts,
                                                             nclustrows,
                                                             nclustcols,
                                                             None)
    if cluster_locs is None:
        cluster_locs = [0] * nclusts
    if row_locs is None:
        row_locs = [0] * nclusts
    if col_locs is None:
        col_locs = [0] * nclusts

    if cluster_scales is None:
        cluster_scales = [default_scale] * nclusts
    if row_scales is None:
        row_scales = [default_scale] * nclusts
    if col_scales is None:
        col_scales = [default_scale] * nclusts
    if error_scales is None:
        error_scales = [0] * nclusts

    row_matrix = _make_row_matrix_(nrows, nclusts, nclustrows, noverlap_rows)
    col_matrix = _make_col_matrix_(ncols, nclusts, nclustcols, noverlap_cols)

    data = numpy.zeros((nrows, ncols)) + dist(loc=background_loc,
                                              scale=background_scale)
    for rows, cols, clust_loc, clust_scale, \
            row_loc, row_scale, col_loc, col_scale \
            in zip(row_matrix.T,
                col_matrix.T,
                cluster_locs,
                cluster_scales,
                row_locs,
                row_scales,
                col_locs,
                col_scales):
        clust_effect = dist(loc=clust_loc, scale=clust_scale)
        row_effects = dist(loc=row_loc, scale=row_scale, size=nrows)
        col_effects = dist(loc=col_loc, scale=col_scale, size=ncols)
        cluster = numpy.vstack(row_effects) + col_effects + clust_effect
        cluster[numpy.outer(rows, cols) <= 0] = 0
        data +=  cluster

    if noise > 0:
        data = add_noise(data, scale=noise)
    data = _add_bicluster_noise_(data, row_matrix, col_matrix, error_scales)

    expected = _make_expected_biclusters_(row_matrix, col_matrix, data)
    if shuffle:
        data, expected = _shuffle_(data, expected)
    return data, expected


def make_fabia_data(nrows,
                    ncols,
                    nclusts,
                    f1,
                    f2,
                    of1,
                    of2,
                    sd_noise,
                    sd_z_noise,
                    mean_z,
                    sd_z,
                    sd_l_noise,
                    mean_l,
                    sd_l,
                    shuffle=True,
                    pos=False):
    """
    Make FABIA-style data.

    An interface to the Bioconductor 'fabia' library's
    makeFabiaDataset functions.

    Requires that 'fabia' be installed.

    Args:
        * nrows: number of observations.
        * ncols: number of samples.
        * nclusts: number of biclusters.
        * f1: ncols/f1 max. additional samples are active in a bicluster.
        * f2: nrows/f2 max. additional observations that form a pattern
            in a bicluster.
        * of1: minimal active samples in a bicluster.
        * of2: minimal observations that form a pattern in a bicluster.
        * sd_noise: Gaussian zero mean noise std on data matrix.
        * sd_z_noise: Gaussian zero mean noise std for deactivated hidden factors.
        * mean_z: Gaussian mean for activated factors.
        * sd_z: Gaussian std for activated factors.
        * sd_l_noise: Gaussian zero mean noise std if no observation patterns
            are present.
        * mean_l: Gaussian mean for observation patterns.
        * sd_l: Gaussian std for observation patterns.
        * shuffle: If True, shuffle dataset.
        * pos: Use the MakeFabiaDataPos functions

    """
    robjects.r.library('fabia')

    function = 'makeFabiaData'
    if not shuffle:
        function += "Blocks"
    if pos:
        function += "Pos"
    func = robjects.r[function]

    result = func(nrows,
                  ncols,
                  nclusts,
                  f1,
                  f2,
                  of1,
                  of2,
                  sd_noise,
                  sd_z_noise,
                  mean_z,
                  sd_z,
                  sd_l_noise,
                  mean_l,
                  sd_l)

    noisy_data = numpy.array(result[0]).copy()
    noiseless_data = numpy.array(result[1]).copy()
    cols_vector = result[2]
    rows_vector = result[3]

    f = lambda x: int(x) - 1
    rows = []
    for r in rows_vector:
        rows.append(map(f, r))

    cols = []
    for c in cols_vector:
        cols.append(map(f, c))

    biclusters = []
    for r, c in zip(rows, cols):
        biclusters.append(Bicluster(r, c, noisy_data))

    return noisy_data, biclusters


def make_isa_data(nrows=300,
                  ncols=50,
                  nclusts=3,
                  nclustrows=None,
                  nclustcols=None,
                  noise=0,
                  bicluster_signals=None,
                  bicluster_noise=None,
                  noverlap_rows=0,
                  noverlap_cols=None,
                  shuffle=None):
    """
    Make ISA-style data.

    Generates a dataset using the Bioconductor 'isa2' package's
    make.isa.data function.

    If an argument is None, it is not included, and isa2's defaults are used.

    Requires that 'isa2' be installed.

    Args:
        * nrows: Number of rows in the data matrix.
        * cols: Number of columns in the data matrix.
        * nclusts: Number of biclusters.
        * nclustrows: Rows in each bicluster.
            Defaults to round(0.5 * num_rows/num_fact)
        * nclustcols: Cols in each bicluster. round(0.5 * num_cols/num_fact)
        * noise: Standard deviation of normal noise in background.
        * bicluster_signals: List of base signals for each bicluster.
            Defaults to 1's.
        * bicluster_noise: List of noise standard deviations for each bicluster.
            Defaults to 0's.
        * noverlap_rows: Number of bicluster rows that overlap.
        * noverlap_cols: Number of coluster columns that overlap.
            Defaults to 'overlap_row'.
        * shuffle: If True, shuffle rows and columns.

    """
    args = locals()

    isa_map = dict(
        nrows='num_rows',
        ncols='num_cols',
        nclusts='num_fact',
        nclustrows='mod_row_size',
        nclustcols='mod_col_size',
        noise='noise',
        bicluster_signals='mod_signal',
        bicluster_noise='mod_noise',
        noverlap_rows='overlap_row',
        noverlap_cols='overlap_col',
        )

    isa_args = dict()

    for key, argkey in isa_map.iteritems():
        isa_args[argkey] = args[key]

    #remove empty keys
    empty_keys = []
    for key in isa_args:
        if isa_args[key] is None:
            empty_keys.append(key)
    for key in empty_keys:
        isa_args.pop(key)

    for key in ['mod_signal', 'mod_noise']:
        if key in isa_args:
            isa_args[key] = robjects.FloatVector(list(isa_args[key]))

    robjects.r.library('isa2')

    #get data
    func = robjects.r['isa.in.silico']
    result = func(**isa_args)

    #convert to python
    data = numpy.array(robjects.Matrix(result[0])).copy()
    rows = numpy.array(robjects.Matrix(result[1])).copy()
    cols = numpy.array(robjects.Matrix(result[2])).copy()

    nbiclusters = rows.shape[1]

    row_list = []
    for i in range(nbiclusters):
        row = list(rows[:,i].nonzero()[0])
        row_list.append(row)

    col_list = []
    for i in range(nbiclusters):
        col = list(cols[:,i].nonzero()[0])
        col_list.append(col)

    expected = []
    for r, c, in zip(row_list, col_list):
        expected.append(Bicluster(r, c, data))

    if shuffle:
        data, expected = _shuffle_(data, expected)
    return data, expected
