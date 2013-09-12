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

"""Utilities for reading and writing datasets for various algorithms."""

from __future__ import division
import numpy as np

def write_expression_data(data,
                          filename,
                          sep='\t',
                          genes=None,
                          conditions=None):
    """Writes a dataset in the following relatively standard format::

        Genes/Conditions [col ID] [col ID] ... [col ID]
        [row ID] [value] [value] ... [value]
        [row ID] [value] [value] ... [value]
        ...
        [row ID] [value] [value] ... [value]

    Args:
        * data: numpy.ndarray
        * filename: Output file name.
        * sep: Seperating character, e.g. ' ' or ','.
        * genes: Optional list of row labels.
        * conditions: Optional list of column labels.

    """
    outfile = file(filename, 'w')
    nrows, ncols = data.shape

    if genes is None:
        genes = ['row{0}'.format(i) for i in range(nrows)]
    else:
        assert len(genes) == nrows

    if conditions is None:
        conditions = ['col{0}'.format(i) for i in range(ncols)]
    else:
        assert len(conditions) == ncols

    outfile.write("Genes/Conditions")
    outfile.write(sep)
    outfile.write("\t".join(conditions))
    for i, gene in enumerate(genes):
        outfile.write("\n")
        outfile.write(gene)
        outfile.write(sep)
        outfile.write(sep.join([str(j) for j in data[i]]))
    outfile.close()


class ExpressionArray(np.ndarray):
    """
    A numpy array with extra attributes 'genes', 'samples', and 'annotation'.

    adapted from http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

    """
    def __new__(cls, input_array, genes=None, samples=None, annotation=None):
        obj = np.asarray(input_array).view(cls)
        obj.genes = genes
        obj.samples=samples
        obj.annotation=annotation
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.genes = getattr(obj, 'genes', None)
        self.samples = getattr(obj, 'samples', None)
        self.annotation = getattr(obj, 'annotation', None)

    def __reduce__(self):
        object_state = list(np.ndarray.__reduce__(self))
        subclass_state = (self.genes, self.samples, self.annotation)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        nd_state, own_state = state
        np.ndarray.__setstate__(self, nd_state)
        self.genes, self.samples, self.annotation = own_state



#TODO: refactor to use labelled arrays or numpy's DataArrays
def read_expression_data(filename, skip_header=1, strip_chars=None):
    """
    Read a tsv file with the same format written by write_expression_data().

    Args:
        * filename:
        * skip_header:
        * strip_chars:

    Returns:
        An instance of ExpressionArray.

    """
    genes = list(np.genfromtxt(filename,
                               skip_header=skip_header,
                               usecols=[0],
                               dtype=np.character))
    data = np.genfromtxt(filename, skip_header=skip_header)[:,1:]
    headers = []
    if skip_header > 0:
        headers = open(filename).read().split('\n')[0:skip_header]
        headers = [line.split() for line in headers]
    if strip_chars is not None:
        genes = [g.strip(strip_chars) for g in genes]
        headers = [[h.strip(strip_chars) for h in header] for header in headers]
    if len(headers[0]) == data.shape[1] + 1:
        headers = [h[1:] for h in headers] #remove top-left header, if present.
    return ExpressionArray(data, genes, headers)


def write_bicoverlapper(bicluster_sets, filename, rownames=None, colnames=None):
    """
    Writes every bicluster in a list of lists of biclusters to a file
    in the format read by BicOverlapper.

    Args:
        * bicluster_sets: a list of list of biclusters that will be
            written to the file.
        * filename: output file name.

    File format::

        [number_of_biclusters]
        bicluster set 1
        #rows bic1.1 #columns bic1.1
        row1 row2 ... rowN
        col1 col2 ... colN
        #rows bic1.2 #columns bic1.2
        row1 row2 ... rowN
        col1 col2 ... colN
        ...
        bicluster set 2
        #rows bic2.1 #columns bic2.1
        row1 row2 ... rowN
        col1 col2 ... colN
        #rows bic2.2 #columns bic2.2
        row1 row2 ... rowN
        col1 col2 ... colN
        ...

    """
    data = bicluster_sets[0][0].data
    for biclusters in bicluster_sets:
        for b in biclusters:
            if not id(b.data) == id(data):
                raise Exception('Bicluster data is not identical.')

    if rownames is None:
        rownames = ['row{0}'.format(i) for i in range(data.shape[0])]
    if colnames is None:
        colames = ['col{0}'.format(j) for j in range(data.shape[1])]
    outfile = file(filename, 'w')
    outfile.write("{0}\n".format(sum([len(b) for b in bicluster_sets])))
    for i, biclusters in enumerate(bicluster_sets):
        outfile.write("bicluster_set_{0}\n".format(i))
        for j, bicluster in enumerate(biclusters):
            outfile.write("#rows_bic{0}.{1} #columns_bic{0}.{1}\n".format(i, j))
            outfile.write(" ".join([rownames[r] for r in bicluster.rows]))
            outfile.write("\n")
            outfile.write(" ".join([colnames[c] for c in bicluster.cols]))
            outfile.write("\n")
    outfile.close()


def write_pcl_dataset(data, filename):
    """
    Given the pure numpy data matrix with only expression values,
    converts the data matrix into PCL format with default row and column names.

    Args:
        * data: numpy.ndarray
        * filename: output file name.

    """
    nrows, ncols = data.shape
    line1 = ["ID", "NAME", "GWEIGHT"]
    line1.extend(range(ncols))

    line2 = ["EWEIGHT", "", "",]
    line2.extend([1] * ncols)

    pcllist = [line1, line2]
    for i, line in enumerate(data):
        pclline = [i, '', 1]
        pclline.extend(line)
        pcllist.append(pclline)

    pcl = '            \n'.join(['\t'.join(map(str, line)) for line in pcllist])
    with open(filename, 'w') as f:
        f.write(pcl)


def write_david_multilist(filename, gene_lists, name=None):
    """
    Writes a DAVID multilist, with each list in one column.

    The first row gives the name of the list, which is just 'name#'.

    The gene names must be the same for each bicluster

    """
    if name is None:
        name = "unnamed_list"
    longest = max([len(i) for i in gene_lists])

    def extend(gene_list, length):
        gl = gene_list[:]
        gl.extend(['' for i in range(length - len(gl))])
        return gl

    gene_lists = [extend(g, longest) for g in gene_lists]
    lines = [[gene_lists[i][j]
              for i in range(len(gene_lists))]
             for j in range(longest)]
    with open(filename, 'w') as f:
        f.write('\t'.join([name + str(i) for i in range(len(gene_lists))]))
        f.write('\n')
        for line in lines:
            f.write('\t'.join(line))
            f.write('\n')


def write_david_list(filename, gene_list):
    """
    Write a DAVID (http://david.abcc.ncifcrf.gov/) list of genes.

    """
    with open(filename, 'w') as f:
        f.write('\n'.join(gene_list))
