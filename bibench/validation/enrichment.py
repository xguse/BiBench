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
Gene Ontology enrichment analysis. Computes the hypergeometric
probability of a set of genes being enriched with a GO term.

This module allows for multiple test correction, but otherwise makes a
lot of assumptions and is rather fragile.  For a more customized
workflow, you will have to write your own enrichment workflow.

Since this module depends on GOStats, it requires that a '.db'
annotation package be available for the organism you are testing. If
no such annotation package is available, try topGO instead.

TODO:
    * Extend to more than just gene ontoloy, e.g. KEGG.
    * Expose more of the R packages' functionality.

"""

import rpy2
from rpy2 import robjects
#from rpy2.interactive import importr # WAD: this does not seem to be the correct way to get to importr anymore
from rpy2.robjects.packages import importr # WAD: this works for me
from rpy2.rinterface import RRuntimeError
from collections import namedtuple
from bibench.rutil import get_bioclite

ontologies = ["BP", "CC", "MF"]
correction_methods = ["Bonferroni",
                      "Holm",
                      "Hochberg",
                      "SidakSS",
                      "SidakSD",
                      "BH",
                      "BY",
                      "ABH",
                      "TSBH"]

EnrichedGO = namedtuple('EnrichedGO', 'goid raw_p corrected_p')

def enrichment(bicluster,
               annotation_db_name,
               gene_universe,
               plimit=0.05,
               ontology='BP',
               mtc='BH',
               conditional=True,
               test_direction='over'):
    """Computes the gene enrichment of the rows of a bicluster.

    Args:
        * bicluster: A Bicluster instance; its rows are tested for enrichment.
        * annotation_db_name: string; the name of an annotation package supported
            by GOStats. e.g.: 'hgu95av2.db'
        * gene_universe: A list of gene names; the gene universe out of which the
            'genes' list was chosen.
        * plimit: The threshold p-value for significance.
        * ontology: Which ontology to use: one of 'BP', 'CC', 'MF'.
        * mtc: Multiple test correction method to use. Options:
            * 'Bonferroni'
            * 'Holm'
            * 'Hochberg'
            * 'SidakSS' (Sidak single-step)
            * 'SidakSD' (Sidak step-down)
            * 'BH' (Benjamini & Hochberg)
            * 'BY' (Benjamini & Yekutieli)
            * 'ABH' (adaptive Benjamini & Hochberg)
            * 'TSBH' (two-stage Benjamini & Hochberg)
        * conditional: If True, use the GO structure. If false, perform a
            standard hypergeometric test.
        * test_direction: 'over' or 'under'. Whether to detect over or under
            enriched terms.

    Returns:
        A list of EnrichedGO tuples, one for each gene ontology
        term that is enriched in this set of genes.

    """
    genes = [gene_universe[i] for i in bicluster.rows]
    return gene_enrichment(genes,
                           annotation_db_name,
                           gene_universe,
                           plimit,
                           ontology,
                           mtc)


def gene_enrichment(genes,
                    annotation_db_name,
                    gene_universe,
                    plimit=0.05,
                    ontology='BP',
                    mtc='BH',
                    conditional=True,
                    test_direction='over'):
    """
    Finds the p-value and enriched GO terms for a set of gene names.

    Args: (see enrichment())
        * genes: A list of gene names to test.

    Returns:
        A list of EnrichedGO tuples, one for each gene ontology
        term that is enriched in this set of genes.

    """
    if not ontology in ontologies:
        raise Exception(
            'Unknown ontology: {0}.' \
                ' Valid arguments: {1}'.format(ontology,
                                               ", ".join(ontologies)))
    if not mtc in correction_methods:
        raise Exception(
            'Unknown multiple test correction method: {0}.' \
                'Valid arguments: {1}'.format(mtc,
                                              ", ".join(correction_methods)))
    if not set(genes).issubset(set(gene_universe)):
        raise Exception("'genes' is not a subset of 'gene_universe'")


    importr('multtest')
    importr('MASS')
    importr('GOstats')

    #load correct annotation library
    if annotation_db_name[-3:] != '.db':
        annotation_db_name = annotation_db_name + '.db'
    try:
        robjects.r.library(annotation_db_name)
    except RRuntimeError:
        bioclite = get_bioclite()
        bioclite(annotation_db_name)
        robjects.r.library(annotation_db_name)

    #convert gene names to R data structures
    genes = robjects.StrVector(genes)
    gene_universe = robjects.StrVector(gene_universe)

    entrez_id = robjects.r.get(annotation_db_name[:-3] + 'ENTREZID')

    selected_entrez_ids = robjects.r.unlist(robjects.r.mget(genes,
                                                            entrez_id,
                                                            ifnotfound=robjects.NA_Character))
    entrez_universe = robjects.r.unlist(robjects.r.mget(gene_universe,
                                                        entrez_id,
                                                        ifnotfound=robjects.NA_Character))

    GOparams = robjects.r.new('GOHyperGParams',
                              geneIds=selected_entrez_ids,
                              universeGeneIds=entrez_universe,
                              annotation=annotation_db_name,
                              ontology=ontology,
                              pvalueCutoff=plimit,
                              conditional=conditional,
                              testDirection=test_direction)

    hg_over = robjects.r['hyperGTest'](GOparams)
    p_over_categories = robjects.r['sigCategories'](hg_over, 1)
    length = len(p_over_categories)
    if length == 0:
        return robjects.DataFrame({})

    #multiple test correction
    pvalues = rpy2.rinterface.globalenv.get("pvalues")
    mt = robjects.r['mt.rawp2adjp'](pvalues(hg_over)[0:length], proc= mtc)

    result =  [EnrichedGO(gs, p, pc)
               for gs, p, pc in zip(p_over_categories,
                                    pvalues(hg_over)[0:length],
                                    mt.rx('adjp')[0].rx(True, 2))]

    #TODO: shouldn't plimit already filter results?
    return [r for r in result if r.corrected_p <= plimit]


GoAnnot = namedtuple('GoAnnot', 'goid term ontology synonym secondary definition')

def goid_annot(goid):
    """
    Get GO annotation for a GO id. If the annotation is empty, you probably need
    to update your GO.db package.

    Depends on the Bioconductor GO.db package.

    >>> get_go_terms('GO:0051649')
    result

    """
    importr('GO.db')
    term = str(robjects.r['Term'](goid)[0])
    ont = str(robjects.r['Ontology'](goid)[0])

    syn = robjects.r['Synonym'](goid)[0]
    if robjects.r['is.null'](syn):
        syn = None
    else:
        syn = [str(i) for i in syn]

    sec = robjects.r['Secondary'](goid)[0]
    if robjects.r['is.null'](sec):
        sec = None
    else:
        sec = [str(i) for i in sec]

    defn = str(robjects.r['Definition'](goid)[0])
    return GoAnnot(goid, term, ont, syn, sec, defn)
