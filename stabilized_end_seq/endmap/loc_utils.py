import numpy as np
from rend_utilities.peak_parsing_utils import is_in_coding, is_near_coding
from rend_utilities.wig_utils import gene_from_position
import pdb

def count_cleavage_sites(start, end, strand, cs_df):
    peaks = cs_df[(cs_df['strand']==strand) & (cs_df['position']>start) & (cs_df['position']<end)]
    return len(peaks)

def get_d_feature_ends(pos, strand, df, strandcol='STRAND', startcol = 'START', endcol = 'END'):
    '''
    Calculate distance of a position to the start and end of a particular feature 
    p: position of peak to consider
    strand: strand of peak
    df: Dataframe encoding gene annotation to check peaks against. Columns from MochiView location set file
    '''
    pos = int(pos)

    # print(strand)
    # if strand == '-':
    #     pdb.set_trace()
    roi = df[(df[strandcol]==strand) & (df[endcol]>pos) & (df[startcol]<pos)]
    if len(roi)>1: # If overlapping genes, multiple ROI will come up. Skip these.
        print('%s on strand %s is present in two overlapping feature windows' % (str(pos),strand))
        return np.nan, np.nan

    if len(roi) == 0: # If peak doesn't fall within a coding region, return nan
        return np.nan, np.nan

    if strand == '+':
        return pos-int(roi[startcol]), int(roi[endcol])-pos
    if strand == '-':
        return int(roi[endcol])-pos, pos-int(roi[startcol])

def get_near_coding(df, plus_genes, minus_genes, search_threshold):
    # Take a dataframe of peaks and add a column specifying whether that peak is near a coding region

    for strand in ('+', '-'):
        df.loc[(df['strand'] == strand), 'near_coding'] = \
            [is_near_coding(pos, strand, search_threshold, plus_genes, minus_genes) for pos in df.loc[df['strand'] == strand,'position']]
    return df

def get_nearby_gene(df, plus_genes, minus_genes):
    df['nearby_gene'] = [None]*len(df)
    for strand in ('+', '-'):
        df.loc[(df['strand'] == strand), 'nearby_gene'] = \
            [str(gene_from_position(pos, strand, plus_genes, minus_genes)) for pos in df[df['strand'] == strand]['position']]
    return df
