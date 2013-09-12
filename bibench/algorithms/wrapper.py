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
Wrapper architecture for tasks common to many biclustering
wrappers. Contains utility functions for writing datasets,
reading back results, etc.

This module is only useful for wrapping binaries that expect data
input as a file, and/or write their results out to a file.

To easily write an interface for such binaries, simply implement
write_dataset(), read_results, and do_call(). Then you can call
wrapper_helper(), passing those functions

To see how this works, try looking at qubic.py, cpb.py, or any other
module that uses this function.

"""

import os
import os.path
import shutil
import tempfile
import subprocess
import logging

import bibench.util as util

class WrapperException(Exception):
    pass


def wrapper_helper(binary,
                   write_dataset,
                   read_results,
                   do_call,
                   data,
                   *args,
                   **kwargs):
    """
    Perform biclustering. Convenience function for wrappers.

    Takes care of creating a temporary directory, exporting the
    dataset, running, and reading back the results.

    Many biclustering implementations require specific data formats
    and produce results in specialized formats. Call function for such
    algorithms, implementing write_dataset(), read_results(), and
    do_call().

    Args:
        * binary: the name of the binary to call.

        * do_call: This function formats the command and runs it. It
          must have the following arguments::

                (data, datafile, results_dir, **kwargs)

          It must call ``binary``, passing it the data in ``datafile``
          and storing results in the ``resultsdir`` directory. Any
          other arguments specific to each algorithm are passed in
          kwargs.

        * write_dataset: The function that writes the dataset out to
          ``filename`` in the necessary format. It must have the
          following arguments::

                (data, filename)

        * read_results: The function that reads the results back and
          returns them as a list of Biclusters. It must have the
          following arguments::

                 (results_dir, data)

     Returns:
         * A BiclusterList of Biclusters.

    """
    if util.which(binary) is None:
        raise WrapperException(
            "executable '{0}' not found on $PATH".format(binary))


    #get a temporary directory to hold dataset and results
    directory = tempfile.mkdtemp()

    #write the dataset in the appropriate format
    datafile = os.path.join(directory, "data.txt")
    write_dataset(data, datafile)

    #prepare location to hold results
    results_dir = os.path.join(directory, "results")
    os.mkdir(results_dir)

    do_call(data, datafile, results_dir, *args, **kwargs)

    #read back the results
    biclusters = read_results(results_dir, data)

    #cleanup the temporary directory
    shutil.rmtree(directory)

    return biclusters
