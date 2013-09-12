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
Classes for represting biclusters, and some utility functions for dealing with
common bicluster tasks, like IO.

"""
from __future__ import division

import copy
from bibench import util
import numpy as np
import inspect
from decorator import decorator


def _get_data_(data1, data2):
    if id(data1) == id(data2):
        return data1
    return None


    biclusters.alg=bbc
    biclusters.args = kwargs


def _merge_dicts_(a, b):
    return dict(a.items() + b.items())


def _get_args_dict_(f, args, kwargs):
    argspec = inspect.getargspec(f)
    firstdefault = -len(argspec.defaults)

    required = argspec.args[:firstdefault]
    optional = argspec.args[firstdefault:]

    req_dict = dict(zip(required, args[:firstdefault]))
    opt_dict = dict(zip(optional, argspec.defaults))

    given_opt_dict = dict(zip(optional[:len(args) - len(required)], args[firstdefault:]))
    final_opt_dict = dict(opt_dict.items() + given_opt_dict.items())

    return dict(req_dict.items() + final_opt_dict.items())


@decorator
def bicluster_algorithm(f, *args, **kwargs):
    """
    Decorator to automatically set 'alg' and 'args' attribute of
    results of a biclustering algorithm.

    """
    result = f(*args, **kwargs)
    props = None
    if hasattr(result, 'properties'):
        props = result.properties
    args_dict = _get_args_dict_(f, args, kwargs)
    fname = '.'.join([f.__module__, f.__name__])
    biclusters = BiclusterList(result, fname, args_dict, props)
    return biclusters


class BiclusterList(list):
    """
    A list of biclusters with three extra attributes:

    * alg: the algorithm that generated these biclusters
    * args: the arguments to 'alg'
    * properties: properties, such as likelihood, of this clustering,
      if any.

    """
    def __init__(self, itr, algorithm=None, arguments=None, properties=None):
        list.__init__(self,itr)
        self.algorithm = algorithm
        self.arguments = arguments
        self.properties = properties


class Bicluster:
    """A class for representing biclusters."""


    def __init__(self, rows, cols, data=None):
        """
        Args:
            * rows: A list of ints; the row indices that make up this bicluster.
            * cols: A list of ints; column indices that make up this bicluster.
            * data: An numpy.ndarray. Dataset on which this bicluster is defined.
                Required by some methods.


        Returns:
            A Bicluster instance.

        """
        self.rows = rows
        self.cols = cols
        self.data = data


    def __eq__(self, other):
        """
        Test two biclusters for equality.

        Two biclusters are equal if the have the same rows and
        columns, and they have the same object as their data
        member. It is not enough that their data be equal; it must be
        the same object.

        Args:
            * other: A bicluster to compare.

        """

        return set(self.rows) == set(other.rows) and \
               set(self.cols) == set(other.cols) and \
               id(self.data) == id(other.data)


    def copy(self):
        """Returns a deep copy of this instance."""
        other = Bicluster(copy.copy(self.rows), copy.copy(self.cols), self.data)
        return other


    def array(self, rows=None, cols=None):
        """
        Get a numpy array bicluster from data, using the indices in bic_indices.

        Note: requires that this Bicluster's data member is not None.

        Args:
            * rows: the row indices to use; defaults to this bicluster's rows.
            * cols: the column indices; defaults to this bicuster's columns.

        """
        if not self.data is None:
            if rows is None:
                rows = self.rows
            if cols is None:
                cols = self.cols
            array = self.data.take(rows, axis=0).take(list(cols), axis=1)
            return array


    def filter_rows(self):
        """
        Returns the dataset with only the rows from this bicluster.

        Note: requires that this Bicluster's data member is not None.

        """
        return self.array(cols=np.arange(self.data.shape[1]))


    def filter_cols(self):
        """
        Returns the dataset with only the columns from this bicluster.

        Note: requires that this Bicluster's data member is not None.

        """
        return self.array(rows=np.arange(self.data.shape[0]))


    def intersection(self, other):
        """
        Returns a new bicluster with common rows and columns.

        Args:
            * other: a Bicluster

        Returns:
            A Bicluster instance, with rows and columns common to both
            self and other.

            If other and self have the same data attribute, the
            returned Bicluster also has it; else its data attribute is
            None.

        """

        rows = set(self.rows).intersection(set(other.rows))
        cols = set(self.cols).intersection(set(other.cols))
        return Bicluster(rows, cols, _get_data_(self.data, other.data))


    def union(self, other):
        """
        Returns a new bicluster with union of rows and columns.

        Args:
            * other: a Bicluster

        Returns:
            A Bicluster instance with all rows and columns from both self
            and other.

            If other and self have the same data attribute, the
            returned Bicluster also has it; else its data attribute is
            None.

        """

        rows = set(self.rows).union(set(other.rows))
        cols = set(self.cols).union(set(other.cols))
        return Bicluster(rows, cols, _get_data_(self.data, other.data))


    def symmetric_difference(self, other):
        """
        Returns a new bicluster with only unique rows and columns,
        i.e. the inverse of the intersection.

        Args:
            * other: a Bicluster

        Returns:
            A Bicluster instance with all rows and columns unique to either self
            or other.

            If other and self have the same data attribute, the
            returned Bicluster also has it; else its data attribute is
            None.

        """

        rows = set(self.rows).symmetric_difference(set(other.rows))
        cols = set(self.cols).symmetric_difference(set(other.cols))
        return Bicluster(rows, cols, _get_data_(self.data, other.data))


    def difference(self, other):
        """
        Returns the difference of two biclusters.

        Args:
            * other: a Bicluster

        Returns:
            A Bicluster instance with self's rows and columns, but not other's.

            If other and self have the same data attribute, the
            returned Bicluster also has it; else its data attribute is
            None.


        """

        rows = set(self.rows).difference(set(other.rows))
        cols = set(self.cols).difference(set(other.cols))
        return Bicluster(rows, cols)



    def issubset(self, other):
        """
        Returns True if self's rows and columns are both subsets of
        other's; else False.

        """

        return (set(self.rows).issubset(set(other.rows)) and
                set(self.cols).issubset(set(other.cols)))


    def shape(self):
        """Returns the number of rows and columns in this bicluster."""
        return len(self.rows), len(self.cols)


    def area(self):
        """Returns the number of elements in this bicluster."""
        return len(self.rows) * len(self.cols)


    def overlap(self, other):
        """Returns the ratio of the overlap area to self's total size."""
        return self.intersection(other).area() / self.area()


    def __repr__(self):
        return "Bicluster({0}, {1})".format(repr(self.rows), repr(self.cols))


def filter(biclusters,
           minrows=2,
           mincols=2,
           max_overlap=1.0,
           remove_subsets=True,
           datashape=None):
    """
    Removes duplicates, small biclusters, overlapping biclusters, and
    biclusters that are as large as the dataset from a list.

    Args:
        * biclusters: a list of biclusters to filter.
        * min_rows: the minimum allowed number of rows.
        * min_cols: the minimum allowed number of columns.
        * max_overlap: the maximum allowed % overlap between any two clusters;
            a float between 0 and 1.
        * remove_subsets: filter out biclusters that are subsets of existing
            biclusters.
        * data: use if bicluster.data is None.

    Returns:
        A sublist of the given biclusters.

    """
    assert minrows > 0 and mincols > 0
    assert max_overlap >= 0 and max_overlap <= 1

    if datashape is not None:
        nrows, ncols = datashape
    else:
        assert np.all([b.data is not None for b in biclusters])
        assert len(set([b.data.shape for b in biclusters])) == 1
        nrows, ncols = biclusters[0].data.shape


    #remove empty biclusters and biclusters as big as the dataset
    biclusters = [b for b in biclusters if len(b.rows) > minrows
                  and len(b.cols) > mincols
                  and (len(b.rows) < nrows
                       or len(b.cols) < ncols)]

    biclusters = sorted(biclusters,
                        cmp = lambda x, y: cmp(x.area(), y.area()),
                        reverse=True)

    def keep(b, accepted):
        for a in accepted:
            if b == a or b.overlap(a) > max_overlap or (remove_subsets
                                                        and b.issubset(a)):
                return False
        return True

    accepted = []
    for b in biclusters:
        if keep(b, accepted):
            accepted.append(b)
    return accepted


def write_biclusters(biclusters, filename):
    """
    Write biclusters to an output file.

    Uses the format:

    <rows>
    <cols>

    seperated by empty lines.

    Args:
        * biclusters: the list of biclusters that will be written to the file
        * filename: a string containing the output file name.

    """

    with open(filename, 'w') as outfile:
        for bicluster in biclusters:
            outfile.write(' '.join([str(r) for r in bicluster.rows]))
            outfile.write("\n")
            outfile.write(' '.join([str(c) for c in bicluster.cols]))
            outfile.write("\n")
            outfile.write("\n")
    outfile.close()


def read_biclusters(filename):
    """Reads the bicluster from a file writtin by write_biclusters().

    Args:
        * filename: a string.

    """
    biclusters = []
    with open(filename) as infile:
        for str_rows, str_cols, newline in util.grouper(infile, 3):
            if str_rows is None or str_cols is None or newline is None:
                continue
            str_rows, str_cols = str_rows.split(), str_cols.split()
            if str_rows == [] or str_cols == []:
                continue
            rows = map(int, str_rows)
            cols = map(int, str_cols)
            biclusters.append(Bicluster(rows, cols))
    return biclusters


def get_row_col_matrices(biclusters):
    """
    Returns the row x number and col x number matrices for the given
    set of biclusters.

    Requires that 'data' member be set and equal for all biclusters.

    Args:
        * biclusters: a list of Bicluster instances.

    Returns:
        The tuple (rowmatrix, colmatrix), where rowmatrix has
        dimensions m by len(biclusters) and colmatrix has dimensions n
        by len(biclusters), where the dataset has m rows and n
        columns.

        Element rowmatrix[x, y] is 1 if row x is in bicluster y, else it is zero.
        Element colmatrix[x, y] is 1 if column x is in bicluster y, else zero.

    """
    data = biclusters[0].data
    assert not data is None
    assert(all([id(b.data) == id(data) for b in biclusters]))

    nrows, ncols = data.shape
    nbiclusters = len(biclusters)

    assert max([max(b.rows) for b in biclusters]) < nrows
    assert max([max(b.cols) for b in biclusters]) < ncols

    RowXNumber = np.zeros((nrows, nbiclusters), dtype=np.bool8)
    ColXNumber = np.zeros((ncols, nbiclusters), dtype=np.bool8)

    for bindex, bicluster in enumerate(biclusters):
        for r in bicluster.rows:
            RowXNumber[r, bindex] = True
        for c in bicluster.cols:
            ColXNumber[c, bindex] = True

    return RowXNumber, ColXNumber
