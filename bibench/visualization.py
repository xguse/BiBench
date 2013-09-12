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

"""Tools for visualizing datasetes and biclustering results."""

import os

from rpy2 import robjects as r
from rpy2.robjects.packages import importr

#enables automatic conversion from numpy to R
import rpy2.robjects.numpy2ri

import pkg_resources
v = pkg_resources.get_distribution('rpy2').version
if v[0:3] >= '2.2':
    rpy2.robjects.numpy2ri.activate()

from bibench.bicluster import get_row_col_matrices


grdevices = importr('grDevices')

devices = dict(png=grdevices.png,
               ps=grdevices.postscript,
               pdf=grdevices.pdf,
               jpg=grdevices.jpeg,
               jpeg=grdevices.jpeg
               )

read_csv = r.r['read.csv'] #use to read a csv file as a data frame

def _get_r_biclust_(biclusters):
    """Takes a list of biclusters and returns an instance of the R Biclust class."""
    #TODO: This is a hacky way to get visualizations.
    r.r.library('biclust')
    classfunc = r.r["BiclustResult"]

    RowxNumber, ColxNumber = get_row_col_matrices(biclusters)
    NumberxCol = ColxNumber.T
    number = len(biclusters)

    empty_list = r.r("list()")
    params = empty_list
    info = empty_list

    return classfunc(empty_list, RowxNumber, NumberxCol, number, info)


def _rplot_(functionname, *args, **kwargs):
    """Utility function for calling plotting functions in R."""
    r.r.library('biclust')
    func = r.r[functionname]

    dkwargs = dict()

    for key in ('file', 'width', 'height'):
        if key in kwargs:
            dkwargs[key] = kwargs.pop(key)

    #TODO: remove code duplication
    extension = None
    if 'file' in dkwargs:
        extension = os.path.splitext(dkwargs['file'])[1][1:]
        device = devices[extension]
        device(**dkwargs)
        func(*args, **kwargs)
        grdevices.dev_off()
    else:
#TODO: get X11 working here, to allow width and height to be set.
#        device = grdevices.X11
#        device(**dkwargs)
        func(*args, **kwargs)


def heatmap(data, bicluster=None, local=False, palette=None, **kwargs):
    """
    Plots the dataset as a heatmap. Optionally rearrange for a bicluster, if supplied.

    Args:
        * data: a numpy.ndarray to plot.
        * bicluster: an optional bicluster
        * local: if True, only plot the bicluster's submatrix; otherwise plot the whole dataset,
            but with rows and columns shuffled so that the bicluster is in the top-left corner.
        * palette: The color palette to use. Must be an RPy2 object representing a color palette,
            eg: palette=rpy2.r.r("heat.colors(10)")
        * kwargs: any other arguments that the R function 'plot' accepts. Common ones include:
            file, width, height.

    """
    kwargs["local"] = local

    #hack to keep visualization from vertical mirroring itself
    if bicluster is None:
        data = data[::-1]

    if bicluster is not None:
        nrows, ncols = data.shape
        assert max(bicluster.rows) < nrows
        assert max(bicluster.cols) < ncols
        kwargs["bicResult"] = _get_r_biclust_([bicluster])
        kwargs["number"] = 1
    if palette is not None:
        kwargs["beamercolor"] = True
        kwargs["paleta"] = palette
    _rplot_("drawHeatmap", r.Matrix(data), **kwargs)


def parallel_coordinates(bicluster,
                         plot='cols',
                         compare=True,
                         info=False,
                         ylab="Value",
                         color=1,
                         **kwargs):
    """
    Parallel coordinate plot of a bicluster.

    Args:
        * bicluster: Bicluster to plot.
        * plot: 'rows', 'cols', or 'both'
        * compare: If True, also plots the rest of the rows/columns in a lighter color.
        * info: If True: include an informative title.
        * ylab: y-axis label.
        * color: foreground color; integer.
        * kwargs: any other arguments that the R function 'plot' accepts. Common ones include:
            file, width, height.

    """
    kwargs.update(locals())

    valid_plots = ['cols', 'rows', 'both']
    if not plot in valid_plots:
        raise Exception("Error: 'plot' argument must be one of: {0}".format(" ".join(valid_plots)))

    for k in ('bicluster', 'plot', 'kwargs', 'color'):
        kwargs.pop(k)

    kwargs['col'] = color

    bicResult = _get_r_biclust_([bicluster])
    number = 1

    assert bicluster.data is not None
    data = r.Matrix(bicluster.data)

    if plot == 'rows':
        kwargs["plotcol"] = False
    elif plot == 'both':
        kwargs["plotBoth"] = True

    _rplot_("parallelCoordinates", data, bicResult, number, **kwargs)


def bubbleplot(data,
               biclusters1,
               biclusters2=None,
               biclusters3=None,
               projection='mean',
               show_labels=False,
               **kwargs):
    """
    A bubbleplot comparison of multiple sets of biclusters which attempts to project them
    down to two dimensions.

    Args:
        * data: numpy.ndarray on which all biclusters are defined.
        * biclusters1: a list of biclusters.
        * biclusters2: a list of biclusters.
        * biclusters3: a list of biclusters.
        * projection: projection method; one of 'mean', 'isomds', 'cmdscale'.
        * show_labels: if True, label each bicluster in the plot.
        * kwargs: any other arguments that the R function 'plot' accepts. Common ones include:
            file, width, height.


    """
    valid_projections = ['mean', 'isomds', 'cmdscale']
    if not projection in valid_projections:
        raise Exception("Error: 'projection' argument must be one of: {0}".format(" ".join(valid_projections)))

    kwargs['projection'] = projection
    kwargs['showLabels'] = show_labels

    bicResult1 = _get_r_biclust_(biclusters1)
    kwargs['bicResult1'] = bicResult1

    if biclusters2 is not None:
        bicResult2 = _get_r_biclust_(biclusters2)
        kwargs['bicResult2'] = bicResult2

    if biclusters3 is not None:
        bicResult3 = _get_r_biclust_(biclusters3)
        kwargs['bicResult3'] = bicResult3

    _rplot_("bubbleplot", r.Matrix(data), **kwargs)
