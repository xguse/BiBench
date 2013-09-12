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

"""ISA biclustering algorithm, which is provided in the R package 'isa2'."""

from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm
from bibench.util import isiterable
import numpy
import rpy2.robjects as robjects


@bicluster_algorithm
def isa(data,
        thr_row=None,
        thr_col=None,
        no_seeds=100,
        direction=['updown', 'updown']):
    """
    ISA biclustering algorithm.

    Args:
        * data: numpy.ndarray.
        * thr_row: threshold value for rows.
        * thr_col: threshold value for cols.
        * no_seeds: number of seeds to generate biclusters.
        * direction: either 'up' for upregulated,
            'down' for downregulated, 'updown' for both(default).

    Returns:
        A list of biclusters.

    """


    #load the isa library
    robjects.r.library('isa2')

    #get an R object for the data
    r_data = robjects.Matrix(data)

    def handle_threshold(x):
        if x is None:
            x = robjects.r['seq'](1, 3, by=0.5)
        else:
            if not isiterable(x):
                x = [x]
            x = robjects.FloatVector(list(x))
        return x

    thr_row = handle_threshold(thr_row)
    thr_col = handle_threshold(thr_col)

    direction = robjects.StrVector(direction)

    #run biclustering
    func = robjects.r('isa')
    result = func(r_data, thr_row, thr_col, no_seeds, direction)

    #get rowXnumber array
    row_matrix = numpy.array(robjects.Matrix(result[0]))

    #get numberXcolumn array
    col_matrix = numpy.array(robjects.Matrix(result[1]))

    num_biclusters = row_matrix.shape[1]
    assert num_biclusters == col_matrix.shape[1]

    #make list of biclusters
    biclusters = []
    for i in range(num_biclusters):
        row_vals = row_matrix[:, i]
        col_vals = col_matrix[:, i]

        rows = [index for index, elt in enumerate(row_vals) if elt]
        cols = [index for index, elt in enumerate(col_vals) if elt]

        biclusters.append(Bicluster(rows, cols, data=data))

    return biclusters
