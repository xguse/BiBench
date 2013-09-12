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
Utilities for finding and importing data from R and Bioconductor.

"""

import numpy
import zlib
import rpy2.interactive as r
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
import pkg_resources
v = pkg_resources.get_distribution('rpy2').version
if v[0:3] >= '2.2':
    rpy2.robjects.numpy2ri.activate()
from bibench.datasets.io import ExpressionArray
import os
from bibench.util import get_hidden_dir, zdumps, zloads
import cPickle

def get_r_data(dataset, library=None):
    """
    Load and return an R dataset as a
    bibench.datasets.io.ExpressionArray.

    Warning: this will only work for datasets that can automatically
    converted to a numpy.ndarra.

    Args:
        * dataset: Name of the dataset to load.
        * library: the name of the library containing the dataset. If
            None, assume it has the same name as the dataset.

    Returns: An ExpressionArray.

    """
    r.importr(library)
    r.packages.utils.data(dataset)
    rdata = robjects.r[dataset]
    genes = list(r.packages.base.rownames(rdata))
    samples = list(r.packages.base.colnames(rdata))
    ndata = numpy.array(r.packages.base.as_matrix(rdata))
    return ExpressionArray(ndata, genes, samples)


def get_bioc_data(dataset, library=None):
    """
    Load and return a Bioconductor ExpressionSet dataset as a
    bibench.datasets.io.ExpressionArray.

    Args:
        * dataset: Name of the dataset to load.
        * library: the name of the library containing the dataset. If
            None, assume it has the same name as the dataset.

    Returns: An ExpressionArray.

    """
    if library is None:
        library = dataset
    r.importr(library)
    r.importr('Biobase')
    r.packages.utils.data(dataset)
    eset = robjects.r[dataset]
    ndata = numpy.array(robjects.r['exprs'](eset))
    result = ExpressionArray(ndata,
                             genes = list(robjects.r['featureNames'](eset)),
                             samples = list(robjects.r['sampleNames'](eset)),
                             annotation = str(robjects.r['annotation'](eset)[0]))
    return result


def get_gds_data(gdsname, destdir=None, pkl=True):
    """
    Get a GDS dataset from the Gene Expression Omnibus. Requres the
    'GEOquery' Bioconductor package.

    Args:
        * gdsname: The GDS dataset number or name.

        * destdir: Where to store downloaded datasets.
            Defaults to '$HOME/.bibench/gds'

        * pkl: cache the final result as compressed binary, rather
            than parse all the time.

    """
    if destdir is None:
        destdir = get_hidden_dir('gds')
    if not os.path.isdir(destdir):
        raise Exception('{0} is not a directory'.format(destdir))

    try:
        int(gdsname)
        gdsname = ''.join(['GDS', str(gdsname)])
    except ValueError:
        pass

    pklpath = os.path.join(destdir, gdsname + '.pkl')
    if pkl:
        if os.path.exists(pklpath):
            try:
                return zloads(open(pklpath, 'rb').read())
            except Exception as e:
                print '{0}'.format(e.message)

    r.importr('GEOquery')
    gds = robjects.r['getGEO'](gdsname, destdir=destdir)

    gplname = robjects.r['Meta'](gds).rx2('platform')[0]
    if gplname is not robjects.NA_Character and \
            gplname is not None and \
            type(gplname) == type('') and \
            gplname[0:3] == 'GPL':
        gpl = robjects.r['getGEO'](gplname, destdir=destdir)
    else:
        gpl = robjects.NULL

    eset =  robjects.r['GDS2eSet'](gds, GPL=gpl)
    annotation = get_annotation(gdsname)

    if annotation is robjects.NA_Character:
        annotation = None
    meta = robjects.r['Meta'](gds)
    ndata = numpy.array(robjects.r['exprs'](eset))
    result = ExpressionArray(ndata,
                             genes = list(robjects.r['featureNames'](eset)),
                             samples = list(robjects.r['sampleNames'](eset)),
                             annotation = annotation)
    if pkl:
        with open(pklpath, 'wb') as f:
            f.write(zdumps(result))
    return result


def get_annotation(gdsname):
    return geo_query(
        'select gpl.bioc_package from ' \
            'gds join gpl on gds.gpl=gpl.gpl ' \
            'where gds="{0}"'.format(gdsname)).rx2(1)[0]


def geo_query(sqlquery):
    """
    Query the SQLlite database of GEO metadata, using any SQL query.

    Downloads the database if it is not already present in $HOME/.bibench.

    """
    r.importr('GEOmetadb')
    destdir = get_hidden_dir()
    sqlname = "GEOmetadb.sqlite"
    sqlpath = os.path.join(destdir, sqlname)
    if not os.path.exists(sqlpath):
        robjects.r['getSQLiteFile'](destdir = destdir)
    con =  robjects.r['dbConnect'](robjects.r['SQLite'](), sqlpath)
    return robjects.r['dbGetQuery'](con, sqlquery)
