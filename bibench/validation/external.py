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
Metrics for evaluating biclusters and sets of biclusters.

"""

from __future__ import division

from collections import namedtuple

import numpy as np

class ExternalError(Exception):
    pass

class EmptyList(ExternalError):
    pass

ListScore = namedtuple('ListScore', 'relevance recovery')

##############################################################
#utilities for union/intersection/difference of areas
#
#note: different from area of union/intersection/difference
##############################################################

def _union_area_(bic1, bic2):
    """
    Returns the area of the element-wise union of biclusters.
    Note: this is different from bic1.union(bic2).area()

    Args:
        * bic1: Bicluster instance
        * bic2: Bicluster instance

    """
    area1 = bic1.area()
    area2 = bic2.area()
    intersection_area = bic1.intersection(bic2).area()
    return area1 + area2 - intersection_area


def _intersection_area_(bic1, bic2):
    """
    Returns the area of the element-wise intersection of biclusters.

    Args:
        * bic1: Bicluster instance
        * bic2: Bicluster instance

    """
    return bic1.intersection(bic2).area()


def _difference_area_(bic1, bic2):
    """
    Returns the element-wise area of bic1 - bic2

    Note: this is different from bic1.difference(bic2).area()

    Args:
        * bic1: Bicluster instance
        * bic2: Bicluster instance

    """
    return bic1.area() - bic1.intersection(bic2).area()


def _area_is_nonzero_(bicluster):
    """
    Raise an exception if bicluster's area is zero, because we need to
    divide by it.

    Args:
        * bicluster: Bicluster instance

    """
    if bicluster.area() == 0:
        nrows, ncols = len(bicluster.rows), len(bicluster.cols)
        msg = 'bicluster is ill defined:' \
            '{0} rows and {1} columns'.format(nrows, ncols)
        raise ExternalError(msg)



################################################
### code for comparing individual biclusters ###
################################################

def recovery(expected, found):
    r"""
    The proportion of true positives found.

    .. math:: \frac{|e \cap f|}{|e|}

    Args:
        * expected: Bicluster.
        * found: Bicluster.

    """
    _area_is_nonzero_(expected)
    return _intersection_area_(expected, found) / expected.area()


def relevance(expected, found, nelts=None):
    r"""
    The proportion of true negatives to all negatives.

    .. math:: \frac{n - |e \cup f|}{n - |e|}

    Args:
        * expected: Bicluster
        * found: Bicluster

        * nelts: the number of elements in the data matrix. This
            parameter can be ommitted if ``expected`` and ``found``
            have their .data attribute set.

    """
    if nelts is None:
        assert expected.data is not None #because we need to know negatives
        assert type(expected.data) == np.ndarray
        nelts = expected.data.size

    true_negs = nelts - _union_area_(expected, found)
    all_negs = nelts - expected.area()

    return true_negs / all_negs


def modified_relevance(expected, found):
    r"""
    The proportion of 'expected' retrieved in 'found'.

    Note that the one-way relevance is not true relevance: it is the
    intersection divided by the size of 'found'.

    One way relevance:

    .. math:: \frac{|e \cap f|}{|f|}

    Args:
        * expected: Bicluster
        * found: Bicluster

    """
    _area_is_nonzero_(found)
    return _intersection_area_(expected, found) / found.area()


def f_measure(expected, found, beta=1, modified=True):
    r"""
    Calculates the f_measure score of a result of an algorithm with
    the expected biclusters that are embedded to the dataset. Uses
    relevance and recovery scores, and takes the harmonic mean of
    them. In the harmonic mean formula, relevance score is scaled
    with the square of beta as in the following formula:

    .. math:: (1 + \beta^2) * (spec * sens) / (\beta^2 * spec + sens)

    Args:
        * expected: Bicluster.
        * found: Bicluster.
        * beta: scale factor in the formula of f_measure.

    """
    if modified:
        spec = modified_relevance(expected, found)
    else:
        spec = relevance(expected, found)
    sens = recovery(expected, found)
    if spec == sens == 0: return 0
    return (1 + beta**2) * (spec * sens) / (beta**2 * spec + sens)


def row_jaccard(expected, found):
    """
    The Jaccard coefficient of rows only.
    As used by Prelic - intersection_of_rows / union_of_rows

    Args:
        * expected: Bicluster
        * found: Bicluster

    """
    intsn = expected.intersection(found)
    union = expected.union(found)
    if len(union.rows) == 0:
        return 0
    return len(intsn.rows) / len(union.rows)


def jaccard(expected, found):
    r"""
    Jaccard coefficient of bicluster area.

    .. math:: \frac{|e \cap f|}{|e \cup f|}

    Args:
        * expected: Bicluster
        * found: Bicluster

    """
    return _intersection_area_(expected, found) / _union_area_(expected, found)



#############################################
### code for comparing lists of biclusters ###
#############################################

def _check_list_(blist):
    """
    Checks if any of the expected or found bicluster lists are empty.

    Args:
        * blist: List of biclusters.

    """
    if len(blist) == 0:
        raise EmptyList('list of biclusters is empty')


def _asym_scores_(expected, found, rel_f, rec_f, comp=max):
    r"""
    Calculates both the relevance and the recovery scores of the
    algorithm result. Each bicluster score is computed according to
    the measure function.  The overall score is computed as the mean
    of the best scores for each bicluster.

    .. math::

        Recovery = \sum_{e \in expected} max_{f \in found} rec\_f(e, f)

        Relevance = \sum_{f \in found} max_{e \in expected} rec\_f(e, f)

    Args:
        * expected: List of target biclusters.
        * found: List of recovered biclusters.
        * rel_f: Function for calculating relevance.
        * rec_f: Function for calculating recovery.
        * comp: Either max() or min(), depending on whether rel_f and
            rec_f should be maximised or minimized.

    Returns: (recovery, relevance) tuple.

    """
    _check_list_(expected)
    _check_list_(found)
    return ListScore(recovery=np.mean([comp([rec_f(e, f)
                                             for f in found])
                                       for e in expected]),
                     relevance=np.mean([comp([rel_f(e, f)
                                              for e in expected])
                                        for f in found]))


def _sym_scores_(expected, found, f, comp=max):
    """
    Same as _asym_scores_, but f is symmetric, and so the order of the
    bicluster arguments does not matter.

    Args:
        * f: a function taking two biclusters and returning a score

    """
    return _asym_scores_(expected, found, f, f)


def prelic_list(expected, found):
    """
    Calculates both the relevance and the recovery scores of the
    algorithm result according to Prelics row jaccard score.

    Args:
        * expected: list of target biclusters
        * found: list of biclusters for comparison

    """
    f = row_jaccard
    return _sym_scores_(expected, found, f)


def f_measure_list(expected, found, beta=1, modified=True):
    r"""
    Calculates both the relevance and the recovery scores of the
    algorithm result according to f_measure score with the given beta
    scale.

    .. math::

        Recovery = \sum_{e \in expected} max_{f \in found} fmeasure(e, f)

        Relevance = \sum_{f \in found} max_{e \in expected} fmeasure(e, f)

    Args:
        * expected: list of target biclusters
        * found: list of biclusters for comparison
        * beta: scale factor in the formula of f_measure.

    """
    f = lambda x, y: f_measure(x, y, beta=beta, modified=modified)
    return _sym_scores_(expected, found, f)


def jaccard_list(expected, found):
    r"""
    Recovery and relevance scores of a set of biclusters using Jaccard
    coefficient.

    .. math::

        Recovery = \sum_{e \in expected} max_{f \in found} \frac{|e \cap f|}{|e \cup f|}

        Relevance = \sum_{f \in found} max_{e \in expected} \frac{|e \cap f|}{|e \cup f|}

    Args:
        * expected: list of target biclusters
        * found: list of biclusters for comparison

    """
    f = jaccard
    return _sym_scores_(expected, found, f)


def recovery_relevance_list(expected, found, modified=True):
    r"""
    N.B.: The f-measure calculated as the harmonic mean of the
    recovery and relevance scores returned by this function is
    *different* from the f_measure_list function.

    .. math::

        Recovery = \sum_{e \in expected} max_{f \in found} \frac{|e \cap f|}{|e|}

        Relevance = \sum_{f \in found} max_{e \in expected} \frac{|e \cap f|}{|f|}


    Args:
        * expected: list of target biclusters
        * found: list of biclusters for comparison

    """
    return _asym_scores_(expected, found, modified_relevance, recovery)
