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
Ordered Preserved SubMatrix wrapper.
"""
from bibench.algorithms.wrapper import wrapper_helper
from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm
import os, subprocess
from bibench import util

BINARY = 'opsm.sh'
RESULTFILE = 'opsm_results.txt'

@bicluster_algorithm
def opsm(data, lValue = 10):
    """
    OPSM biclustering algorithm. Finds biclusters that have non-decreasing rows.

    Args:
        * data: 2 dimensional numpy.array format to represent the data matrix.
        * lValue: the number of passed models for each iteration. Default is 10.

    Returns:
        A list of biclusters.

    """
    kwargs = locals()
    return wrapper_helper(BINARY,
                          _write_dataset_,
                          _read_results_,
                          _do_call_,
                          **kwargs)


def _do_call_(data, datafile, results_dir, **kwargs):
    """Executes the OPSM.jar executable with given parameters"""
    outpath = os.path.join(results_dir, RESULTFILE)
    cmd = [BINARY, datafile, str(data.shape[0]),
           str(data.shape[1]), outpath, str(kwargs["lValue"])]
    subprocess.check_call(cmd)


def _read_results_(dirname, data):
    """Read OPSM result file and return list of Bicluster instances."""
    infile = file(os.path.join(dirname, RESULTFILE), 'r')
    geneLine = ""
    conditionLine = ""
    biclusters = []
    for geneLine, conditionLine, seperator in util.grouper(infile, 3):
        bic = _createBicluster_(geneLine, conditionLine, data)
        biclusters.append(bic)
    infile.close()
    return biclusters


def _createBicluster_(geneLine, conditionLine, data):
    """
    Extracts the rows and columns of the bicluster from the given gene
    and condition lines.

    """
    genes = map(int, geneLine.split(" "))
    conditions = map(int, conditionLine.split(" "))
    return Bicluster(genes, conditions, data)


def _write_dataset_(data, filename):
    """
    Writes a dataset in the format for OPSM.
    Separates with tab.
    """
    f = file(filename, 'w')
    for row in data:
        f.write('\t'.join(map(str, row)))
        f.write('\n')
    f.close()
