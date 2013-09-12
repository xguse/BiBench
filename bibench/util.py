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

import os
import itertools
import bibench
import numpy as np
import zlib
import cPickle


isiterable = lambda obj: isinstance(obj, basestring) or getattr(obj, '__iter__', False)

def flatten(nested):
    """
    Flatten a list of lists into a single list.

    >>> flatten([[1, 2, 3], [4, 5, 6]])
    [1, 2, 3, 4, 5, 6]

    """
    return [item for sublist in nested for item in sublist]


def grouper(iterable, n, fillvalue=None):
    """
    Iterate over a list in chunks. From
    'http://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks'

    >>> list(grouper([1, 2, 3, 4], 3, 'x'))
    [(1, 2, 3), (4, 'x', 'x')]

    """
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)


def which(program):
    """
    Check for an executable on the PATH; return its absolute path.

    Taken from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def dict_combinations(d):
    """
    Takes a dictionary containing lists. Generates all combinations of values from those lists.

    Useful for ranges of parameters for functions.

    >>> [i for i in dict_combinations(dict(first=[1,2]))]
    [{'first': 1}, {'first': 2}]

    """
    keys, values = zip(*d.items())
    for x in itertools.product(*values):
        yield dict(zip(keys, x))


def make_index_map(mylist):
    """Map each item in the list to its list index."""
    d = dict()
    for i, item in enumerate(mylist):
        d[item] = i
    return d


def bootstrap(data, size):
    """
    Bootstrap a new dataset, of any size, from the given dataset, with
    replacement.

    Args:
        * data: numpy.ndarray
        * size: int or sequence of ints

    Returns:
        A numpy.ndarray

    """
    choices = np.random.random_integers(low=0, high=data.size-1, size=size)
    return np.array(data.flatten())[choices]


def shuffle(data):
    """Shuffle an array along all axes. Returns the shuffled array."""
    fdata = data.flatten()
    np.random.shuffle(fdata)
    fdata.shape = data.shape
    return fdata

def get_hidden_dir(subdir=None):
    """
    Get the BiBench cache directory, and create it if necessary.

    Args:
        * subdir: a subdir to create if it does not exist.

    """
    destdir = os.path.join(os.getenv("HOME"), '.bibench')
    if not os.path.exists(destdir):
        os.mkdir(destdir)
    if not subdir is None:
        destdir = os.path.join(destdir, subdir)
        if not os.path.exists(destdir):
            os.mkdir(destdir)
    return destdir


def zdumps(obj):
    """dump an object, compressing as much as possible"""
    return zlib.compress(cPickle.dumps(obj,cPickle.HIGHEST_PROTOCOL),9)


def zloads(zstr):
    """load a compressed string dumped by _zdumps_"""
    return cPickle.loads(zlib.decompress(zstr))
