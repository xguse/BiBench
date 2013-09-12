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
Correlated Pattern Bicluster (CPB) algorithm. Finds biclusters with genes
that have large pairwise Pearson correlation.
"""
import os
import subprocess

import bibench
from bibench.algorithms.wrapper import wrapper_helper
from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm, filter

BINARY = 'cpb'
INIT_BINARY = 'init_bicluster'


@bicluster_algorithm
def cpb(data,
        nclus,
        targetpcc=0.9,
        fixed_row=-1,
        fixed_col=-1,
        fixw=0,
        min_seed_rows=3,
        max_seed_rows=None):
    """
    Wrapper for the CPB binary. Finds biclusters with high row-wise correlation.

    Args:
        * data: numpy.ndarray
        * nclus: Number of biclusters to find.
        * targetpcc: Minimum PCC for rows.
        * fixed_row: A row that must be in each bicluster; -1 means none.
        * fixed_col: A column that must be in each bicluster; -1 means none.
        * fixw: Weight for computing error of fixed rows.
        * min_seed_rows: Minimum number of rows in each seed bicluster.
        * max_seed_rows: Maximum number of rows in each seed bicluster.

    Returns:
        A list of biclusters.

    """

    nrows, ncols = data.shape
    if max_seed_rows is None:
        max_seed_rows = nrows

    #check args
    assert nclus > 0
    assert targetpcc <= 1 and targetpcc >= 0
    assert (fixed_row >= 0 and fixed_rows < nrows) or fixed_row == -1
    assert (fixed_col >= 0 and fixed_col < ncols) or fixed_col == -1
    assert fixw <= 1 and fixw >= 0
    assert min_seed_rows >= 1
    assert max_seed_rows >= min_seed_rows and max_seed_rows <= nrows


    kwargs = locals()
    biclusters = wrapper_helper(BINARY,
                                write_dataset=_write_dataset_,
                                read_results=_read_results_,
                                do_call=_do_call_,
                                **kwargs)
    return biclusters


def cpb_filter(biclusters,
               data,
               nclus,
               *args,
               **kwargs):
    """
    Filter out small biclusters found by chance. 'nclus' should be
    large enough to generate a representative sample set.

    Args:
        * biclusters: a list of biclusters found by CPB.
        * data: the dataset they were all run on.
        * nclus: the number of clusters to generate for filtering.
        * args: any parameters, in order, that cpb() takes.
        * kwargs: may be any of the same named parameters as cpb() takes.
            For accurate results, use the same parameters used to
            generate the biclusters to be filtered.

    Returns:
        A sublist of 'biclusters', containing only those
        biclusters that are not likely due to random chance.

    """
    data = bibench.util.shuffle(data)
    kwargs['fixed_row'] = -1
    kwargs['fixed_col'] = -1
    shuffle_results = filter(cpb(data, nclus, *args, **kwargs))
    if len(shuffle_results) == 0:
        return biclusters

    maxarea = max([b.area() for b in shuffle_results])
    return [b for b in biclusters if b.area() > maxarea]


def _make_init_file_(initfile,
                     init_binary,
                     nrows,
                     ncols,
                     nclus,
                     min_seed_rows,
                     max_seed_rows,
                     fixed_row,
                     fixed_col,
                     **kwargs):
    kwargs = locals()
    command = "{init_binary} {nrows} {ncols} {nclus} {min_seed_rows} " \
        "{max_seed_rows} {fixed_row} {fixed_col} {initfile}".format(**kwargs)
    try:
        subprocess.check_call(command.split())
    except OSError:
        raise Exception(
            "Error calling '{0}'. Is it on the PATH?"
            .format(kwargs['init_binary']))


def _do_call_(data, datafile, results_dir, **kwargs):
    datafile = os.path.abspath(datafile)
    results_dir = os.path.abspath(results_dir)

    #initial bicluster file
    directory = os.path.split(results_dir)[0]
    kwargs['initfile'] = os.path.join(directory, 'initfile.txt')
    kwargs['init_binary'] = INIT_BINARY
    _make_init_file_(**kwargs)

    #change to the results directory
    saved_path = os.getcwd()
    try:
        os.chdir(results_dir)
        command = '{0} {1} {initfile} 1 {targetpcc} {fixw}'.format(BINARY, datafile, **kwargs)
        subprocess.check_call(command.split())

    except OSError:
        raise Exception("Error calling 'cpb'. Is it on the PATH?")

    finally:
        os.chdir(saved_path)


def _write_dataset_(data, filename):
    """Writes a dataset to the format that CPB reads:

    [num rows] [num columns]
    [first row's values]
    [second row's values]
    ...
    [last row's values]

    """
    outfile = file(filename, 'w')
    rows, cols = data.shape
    outfile.write("{0} {1}".format(rows, cols))
    for row in data:
        outfile.write("\n")
        outfile.write(" ".join([str(j) for j in row]))
    outfile.write("\n")
    outfile.close()


def _read_results_(results_dir, data):
    #find the results file with extension .out
    files = os.listdir(results_dir)
    outfiles = [f for f in files if os.path.splitext(f)[1] == '.out']

    biclusters = []
    for f in outfiles:
        biclusters.append(_read_result_file_(os.path.join(results_dir, f), data))

    return biclusters


def _read_result_file_(filename, data):
    """
    Reads the bicluster in a single CPB output file.

    The file format is:

    ROWS
    [row index]     [row score]
    [row index]     [row score]
    ...
    [row index]     [row score]
    COLS
    [col index]     [col score]
    [col index]     [col score]
    ...
    [col index]     [col score]

    """
    rows, cols = [], []
    with open (filename, 'r') as f:
        target = rows
        for line in f:
            if line[0] == 'R':
                continue
            elif line[0] == 'C':
                target = cols
                continue
            else:
                target.append(int(line.split()[0]))
        rows.sort()
        cols.sort()
    return Bicluster(rows, cols, data=data)
