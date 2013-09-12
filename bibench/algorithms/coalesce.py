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

"""Coalesce algorithm wrapper for finding biclusters with up and down regulated TF."""
import os, subprocess

from bibench.algorithms.wrapper import wrapper_helper
from bibench.bicluster import \
    Bicluster, BiclusterList, bicluster_algorithm
from bibench.datasets import io
import bibench.util as util

BINARY = 'COALESCE'

@bicluster_algorithm
def coalesce(data,
             geneModuleProbability=0.95,
             conditionPvalueThreshold=0.05,
             conditionZThreshold=0.5,
             normalize=False):
    """
    Wrapper for the COALESCE binary.

    Args:
        * data: numpy.ndarray
        * geneModuleProbability: the probability threshhold for including
            genes in a regulatory module.
        * conditionPvalueThreshold: the P-value threshhold for including
            conditions in a regulatory module.
        * conditionZThreshold: the Z-score threshhold for including
            conditions in a regulatory module.
        * normalize: whether to normalize the data.

    Returns:
        A list of biclusters.

    """
    if normalize is False:
        normalize = 0
    else:
        normalize = 1
    kwargs = locals()
    return wrapper_helper(BINARY,
                          _write_dataset_,
                          _read_results_,
                          _do_call_,
                          **kwargs)


def _do_call_(data, datafile, results_dir, **kwargs):
    """Executes the COALESCE with given parameters"""

    command = "{binary} -i {0}" \
        " -p {geneModuleProbability}" \
        " -c {conditionPvalueThreshold}" \
        " -C {conditionZThreshold}".format(datafile, binary=BINARY, **kwargs)

    if kwargs["normalize"] is not 0:
        command += " -e "

    stndout = os.path.join(results_dir, "bic.out")
    stnderr = os.path.join(results_dir, "debug.out")

    with open(stndout, 'w') as out:
        with open(stnderr, 'w') as err:
            subprocess.check_call(command.split(), stdout=out, stderr=err)


def _read_results_(dirname, data):
    """
    Read the result file which is bic.out and returns the list of
    bicluster objects.

    """
    bicOut = os.path.join(dirname, "bic.out")
    f = open(bicOut,'r')

    biclusters = []
    for clusterLine, geneLine, conditionLine, motifLine in util.grouper(f, 4):
        bic = _createBicluster_(geneLine, conditionLine, data)
        biclusters.append(bic)
    f.close()
    return biclusters


def _createBicluster_(geneLine, conditionLine, data):
    """
    Extracts the rows and columns of the bicluster from the given gene
    and condition line

   """
    genes = map(int, geneLine.split('\t')[1:])
    conditions = map(int, conditionLine.split('\t')[1:])
    return Bicluster(genes, conditions, data)


def _write_dataset_(data, filename):
    """Writes a dataset in the format for Coalesce into pcl format."""
    io.write_pcl_dataset(data, filename)
