import numpy as np
from collections import OrderedDict
from plastid import GenomeArray, GenomicSegment
from scipy.stats.mstats import winsorize
import pandas as pd


def rpkm(count_vector, sequencing_depth):
    k = len(count_vector)/1000.
    r = sum(count_vector)
    m = float(sequencing_depth)
    return r / (k * m)

def rend_RPKM(genomearray, plus_genes, minus_genes, trim=30, winsorization_level=0.05):
    '''
    Calculate RPKM for all genes within
    '''
    print("Calculating RPKM...")  # It takes a few seconds to run, so this is used to track if multiple samples

    RPKM_dict = OrderedDict()

    CDS_depth = 0.000001 * calculate_total_CDS_depth(genomearray, plus_genes,
                                                     minus_genes, trim=trim, winsorization_level=winsorization_level)

    for gene in plus_genes:
        gene_counts = np.array(
            genomearray.get(GenomicSegment("NC_000964.3", plus_genes[gene][0], plus_genes[gene][1], "+")))
        gene_counts = gene_counts[trim:-trim]  # Trim ends of CDS to avoid issues associated w/ end-enrichment
        gene_counts = winsorize(gene_counts, winsorization_level)

        RPKM_dict[gene] = rpkm(gene_counts, CDS_depth)

    for gene in minus_genes:
        gene_counts = np.array(
            genomearray.get(GenomicSegment("NC_000964.3", minus_genes[gene][0], minus_genes[gene][1], "-")))
        gene_counts = gene_counts[trim:-trim]  # Trim ends of CDS to avoid issues associated w/ end-enrichment
        gene_counts = winsorize(gene_counts, winsorization_level)
        RPKM_dict[gene] = rpkm(gene_counts, CDS_depth)

    return RPKM_dict