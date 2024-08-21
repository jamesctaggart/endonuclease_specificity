'''
This module contains functions for generating a DataFrame summarizing end sequencing results
This DataFrame will contain locations of putative cleavage sites given 5'/3' end and Rend-seq data for exo knockouts +/- an endonuclease
'ref' in this module refers to an exo-deficient strain without endo knockouts (e.g. CCB434 = ∆rnjA), with 'ko' refering to exo and endo-deficient (e.g. CCB760 = ∆rny, depletion of rnjA)

Functions in this module are used extensively by script generate_endseq_CleavageSummary.py

The central function in this module is build_endseq_table()
In this function, many helper functions create individual DataFrame objects with position, strand, and a newly calculated column
In the case of multiple knockouts being considered, these dataframes also contain a ko_id column, which allows
metrics like local depth in knockout to be separately calculated for each considered knockout strain

build_endseq_table() merges all of these DataFrames together into one large object which can be filtered
to call putative cleavage sites

This resulting dataframe will be 2 x len(genome) x number_of_knockouts rows long, where 2 is derived from 2 strands.

Functions
---------
get_deadzone(em, w, zthresh, d_thresh, gap, av_thresh, em_type, local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05, deadzone_dist=40, save_z = False, outcol='deadzone')
    Calculates a region around each peak to ignore in 5' end analysis due to technical limitations (fragments expected to be lost in size selection)
    Input data (em) is EndMap object containing 3' end mapped Rend-seq data
calc_ratio(ref_em, ko_em, ko_id, ref_suffix = '_ref', ko_suffix = '_ko', countcol = 'counts')
    Calculate position-wise ratio (WT/KO) at each position between two datasets in EndMap objects
get_all_ratios(ref_df, ko_dfs, ko_ids, countcol = 'counts', ref_suffix = '_ref', ko_suffix = '_ko')
    Runs calc_ratio() for all knockouts for a given reference dataset (ko_dfs is a list of EndMap objects, with matching strings in ko_ids)
get_skip_positions
    Returns a dict (keys = strand, '+' or '-') of lists with positions to be skipped in analysis (below a read count threshold)
calc_depth_step
    Given Rend-seq WT/KO data, calculate the step in Rend-seq density and local depth at each position
get_depth_step_rend_filters
    Runs calc_depth_step() for all knockouts for a given reference dataset
    Concatenates results into one long dataframe
get_sequences_from_position
    Given a position and genome seq, returns the local sequence around it
get_site_kmers
    Given endmap, extracts short sequence (k-mer) centered around each position
background_ratio
    Given an array of counts, calculates the ratio of signal at a position to a window centered around it
    with a gap (in case peaks are broad). This signal to background ratio used as peak calling metric.
calc_background_ratio
    Calculates background ratio (bg_ratio) and depth in background window for all non-skipped positions
get_background_ratios
    Runs calc_background_ratio() for all knockouts for a given reference dataset
get_threeprime_peak_ratio
    Given EndMap objects for 5' end and 3' end sequencing data (reference only, not KO)
    calculate background_ratio in 3' end data with a 1 nt offset
    This offset can also be set to 3 nucleotides, which is necessary for EndoA cleavage sites.
calc_zscore
    Given WT/KO EndMap, calculate local z-score in in each
get_zscores
    Runs calc_zscore() for all knockouts for a given reference dataset
get_groups
    Calculates a column which can be used to optionally group together adjacent positions such that
    broad peaks are not counted as multiple cleavage sites.
build_endseq_table
    Main function which uses the above helpers to generate a large DataFrame that can be filtered to
    extract putative positions of cleavage, e.g. by get_cleavage_positions()
    After constructing DataFrame, calculates 'ratio_bg_ratio' which is ratio of WT/KO background ratios
write_sequences_from_lists
    Used to take list of positions write a file that contains sequences around those positions
write_sequences_from_df
    Used to take DataFrame with 'position' and 'strand' columns and write a file that contains sequences around those positions
write_kplogo_meme_from_df
    Uses write_sequences_from_df to make files compatible for input into kpLogo and meme
write_kplogo_meme_from_lists
    Uses write_sequences_from_lists to make files compatible for input into kpLogo and meme
strand_specific_groupby
    Groups together positions with the same group as calculated by get_groups()
get_cleavage_positions
    Filters a DataFrame built by build_endseq_table to extract putative positions of cleavage
    Filters based on 'ref_bg_ratio', 'ref_depth', 'minus1_3end_bg_ratio' and 'deadzone_5endseq'
    Can be told to use 'minus3_3end_bg_ratio' instead for EndoA site analysis.

'''

import numpy as np
from rend_utilities.wig_utils import calculate_local_depth, calculate_step, check_nearby_peak
import pandas as pd
from scipy.stats.mstats import winsorize
from peak_calling.zscore import zscore
import pdb


# Function for calculating deadzone where we can't find peaks in end sequencing due to size selection
def get_deadzone(em, w, zthresh, d_thresh, gap, av_thresh, em_type, local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05, deadzone_dist=40, save_z = False, outcol='deadzone'):
    '''
    Calculate positions which are near peaks such that end sequencing size selections will make it impossible to capture signal 
    '''
    em.call_peaks_zscore(w, zthresh, d_thresh, gap, av_thresh, local_winsorize=local_winsorize, skip_thresh=skip_thresh, winsor_thresh = winsor_thresh, save_z = save_z)

    f_peaks = em.peaktable.loc[em.peaktable['strand']=='+','position'].to_list()
    r_peaks = em.peaktable.loc[em.peaktable['strand']=='-','position'].to_list()

    df = em.readcounts.reset_index()
    # pdb.set_trace()

    if em_type == '5end':
        df.loc[df['strand']=='+',outcol] = df.loc[df['strand']=='+'].apply(lambda x: check_nearby_peak(x['position'], reference_set=f_peaks, threshold=deadzone_dist, side='left'),axis=1)
        df.loc[df['strand']=='-',outcol] = df.loc[df['strand']=='-'].apply(lambda x: check_nearby_peak(x['position'], reference_set=r_peaks, threshold=deadzone_dist, side='right'),axis=1)
    if em_type == '3end':
        df.loc[df['strand']=='+',outcol] = df.loc[df['strand']=='+'].apply(lambda x: check_nearby_peak(x['position'], reference_set=f_peaks, threshold=deadzone_dist, side='right'),axis=1)
        df.loc[df['strand']=='-',outcol] = df.loc[df['strand']=='-'].apply(lambda x: check_nearby_peak(x['position'], reference_set=r_peaks, threshold=deadzone_dist, side='left'),axis=1)

    df = df[['position','strand',outcol]]

    return df


# region functions to calculate ratio of WT to KO
def calc_ratio(ref_em, ko_em, ko_id, ref_suffix = '_ref', ko_suffix = '_ko', countcol = 'counts'):
    # Reset index and strip away all columns but the ones of interest - position/strand and the quantification column
    # Specify countcol to something like CDS-normalized data if you want, or just use raw 'counts'
    print('Calculating ratio for %s... '% ko_id)
    ref_df = ref_em.readcounts.reset_index()
    ko_df = ko_em.readcounts.reset_index()

    ref_df = ref_df[['position','strand',countcol]]
    ko_df = ko_df[['position','strand',countcol]]

    # Make new dataframe using
    out_df = ref_df.merge(ko_df, on=['position','strand'], how='outer', suffixes=[ref_suffix, ko_suffix])
    out_df['ko_id'] = ko_id

    out_df['5end_ratio'] = out_df[countcol+ref_suffix] / out_df[countcol+ko_suffix]

    return out_df

def get_all_ratios(ref_df, ko_dfs, ko_ids, countcol = 'counts', ref_suffix = '_ref', ko_suffix = '_ko'):
    '''
    For each knockout sample in ko_ids, calculate the ratio between KO and a reference WT
    '''
    print('Calculating all ratios...')
    rat = [calc_ratio(ref_df, ko_dfs[i], ko_ids[i], countcol=countcol,ref_suffix=ref_suffix,ko_suffix=ko_suffix) for i in range(len(ko_ids))]
    print('Concatenating ratio table...')
    return pd.concat(rat)
#endregion

#region functions to process rend-seq data and establish depth (and potentially step) thresholds
def get_skip_positions(ref_5endseq_em, thresh, datacol = 'counts'):
    df = ref_5endseq_em.readcounts

    df.loc[df[datacol]>=thresh,'skip'] = False
    df.loc[df[datacol]<thresh,'skip'] = True

    plus = df.loc[df['strand']=='+','skip'].to_list()
    minus = df.loc[df['strand']=='-','skip'].to_list()
    return {'+':plus, '-':minus}

def calc_depth_step(ref_em, ko_em, ko_id, w, gap, skip_positions = None, winsor_level = 0, ref_suffix = '_ref', ko_suffix = '_ko', countcol = 'counts', mode='up'):

    print('Calculating step/depth for %s... '% ko_id)
    ref_df = ref_em.readcounts.reset_index()
    ko_df = ko_em.readcounts.reset_index()

    # Skip positions is defined as the positions with signal that exceed a certain threshold
    ref_df.loc[ref_df['strand']=='+','skip'] = skip_positions['+']
    ref_df.loc[ref_df['strand']=='-','skip'] = skip_positions['-']
    ko_df.loc[ko_df['strand']=='+','skip'] = skip_positions['+']
    ko_df.loc[ko_df['strand']=='-','skip'] = skip_positions['-']

    ref_df_datasets = {'+': ref_df.loc[ref_df['strand']=='+',countcol].to_list(),'-': ref_df.loc[ref_df['strand']=='-',countcol].to_list()}
    ko_df_datasets = {'+': ko_df.loc[ko_df['strand']=='+',countcol].to_list(),'-': ko_df.loc[ko_df['strand']=='-',countcol].to_list()}

    ref_df.loc[ref_df['skip']==False,'ref_step'] = ref_df.loc[ref_df['skip']==False].apply(lambda x: calculate_step(x['position'], x['strand'], ref_df_datasets[x['strand']], w, gap, winsorization_level=winsor_level),axis=1)
    ko_df.loc[ko_df['skip']==False,'ko_step'] = ko_df.loc[ko_df['skip']==False].apply(lambda x: calculate_step(x['position'], x['strand'], ko_df_datasets[x['strand']], w, gap, winsorization_level=winsor_level),axis=1)

    ref_df.loc[ref_df['skip']==False,'ref_depth'] = ref_df.loc[ref_df['skip']==False].apply(lambda x: calculate_local_depth(x['position'], x['strand'], ref_df_datasets[x['strand']], w, gap, winsorization_level=winsor_level, mode=mode),axis=1)
    ko_df.loc[ko_df['skip']==False,'ko_depth'] = ko_df.loc[ko_df['skip']==False].apply(lambda x: calculate_local_depth(x['position'], x['strand'], ko_df_datasets[x['strand']], w, gap, winsorization_level=winsor_level, mode=mode),axis=1)

    ref_df = ref_df[['position','strand','ref_step', 'ref_depth']]
    ko_df = ko_df[['position','strand','ko_step', 'ko_depth']]

    out_df = ref_df.merge(ko_df, on=['position','strand'], how='outer', suffixes=[ref_suffix, ko_suffix])
    out_df['ko_id'] = ko_id

    return out_df

def get_depth_step_rend_filters(ref_rend_em, ko_rend_ems, ko_ids, w, gap, skip_positions = None, winsor_level = 0, countcol = 'counts', mode = 'up'):
    print('Calculating all step/depth...')
    rend = [calc_depth_step(ref_rend_em, ko_rend_ems[i], ko_ids[i], w, gap, winsor_level = winsor_level, countcol=countcol, mode=mode, skip_positions=skip_positions) for i in range(len(ko_ids))]
    print('Concatenating rend threshold tables...')
    return pd.concat(rend)
# endregion

# region functions to calculate local sequence column
def get_sequences_from_position(p, strand, w, genome_seq):
    # print(p)
    if np.isnan(p):
        return np.nan
    if strand == '+':
        return str(genome_seq[p-(w+1):p+(w-1)])
    elif strand == '-':
        return str(genome_seq[p-w:p+w].reverse_complement())

def get_site_kmers(endmap, half_length, genome_seq, colname = 'local_kmer', skip_positions = None):
    print('Extracting kmers...')
    df = endmap.readcounts.reset_index()

    df.loc[df['strand']=='+','skip'] = skip_positions['+']
    df.loc[df['strand']=='-','skip'] = skip_positions['-']

    df.loc[df['skip']==False, colname] = df.loc[df['skip']==False].apply(lambda x: get_sequences_from_position(x['position'], x['strand'], half_length, genome_seq),axis=1)

    df = df[['position','strand','local_kmer']]
    return df
# endregion

# region Functions to calculate peak scores (beyond just counts)
def background_ratio(p, strand, count_array, w, g, winsorization_level=0.0, mode='both', p_offset = 0, return_depth = False):
    '''
    g is gap between peak and start of counting background. Skip g-1 positions up and downstream of the peak (i.e. g=3 skips 2 nucleotides)
    '''
    assert mode in ('up','down','both')
    p = p-1 # adjust for 1/0-base difference between wig format and python
    assert strand in ('+','-')
    if strand == '+':
        p = p-p_offset # If you want to anchor search farther up/downstream (e.g. for 3' peaks from ndoA sites)
    elif strand == '-':
        p = p+p_offset

    if strand == '+':
        up = count_array[p-w-g:p-g]
        down = count_array[p+g+1:p+w+g+1]
    elif strand == '-':
        up = count_array[p+g+1:p+w+g+1]
        down = count_array[p-w-g:p-g]
    if mode == 'both':
        background = np.concatenate([up, down])
    elif mode == 'up':
        background = up
    elif mode == 'down':
        background = down

    if mode=='up' and p<(w+g) and strand == '+': # added 220615 - fix border cases when quantifying only on one wide - will return an empty array at edges when gap<p<w
        return np.nan
    elif mode=='up' and len(count_array)-p<(w+g) and strand == '-':
        return np.nan
    elif mode=='down' and len(count_array)-p<(w+g) and strand == '+':
        return np.nan
    elif mode=='down' and p<(w+g) and strand == '-':
        return np.nan # new until here

    if return_depth:
        return np.mean(winsorize(np.array(background),winsorization_level))
    else:
        return count_array[p] / np.mean(winsorize(np.array(background),winsorization_level))

def calc_background_ratio(ref_em, ko_em, ko_id, skip_positions=None, countcol = 'counts', w=50, g=3, winsorization_level=0.0, mode='down'):
    print('Calculating peak ratios to background...')

    ref_df = ref_em.readcounts.reset_index()
    ko_df = ko_em.readcounts.reset_index()

    # Skip positions is defined as the positions with signal that exceed a certain threshold
    ref_df.loc[ref_df['strand']=='+','skip'] = skip_positions['+']
    ref_df.loc[ref_df['strand']=='-','skip'] = skip_positions['-']
    ko_df.loc[ko_df['strand']=='+','skip'] = skip_positions['+']
    ko_df.loc[ko_df['strand']=='-','skip'] = skip_positions['-']

    ref_datasets = {'+': ref_df.loc[ref_df['strand']=='+',countcol].to_list(),'-': ref_df.loc[ref_df['strand']=='-',countcol].to_list()}
    ko_datasets = {'+': ko_df.loc[ko_df['strand']=='+',countcol].to_list(),'-': ko_df.loc[ko_df['strand']=='-',countcol].to_list()}

    ref_df.loc[ref_df['skip']==False, 'ref_bg_ratio'] = ref_df.loc[ref_df['skip']==False].apply(lambda x: background_ratio(x['position'],x['strand'],\
        ref_datasets[x['strand']],w=w,g=g,winsorization_level=winsorization_level,mode=mode), axis=1)
    ko_df.loc[ko_df['skip']==False, 'ko_bg_ratio'] = ko_df.loc[ko_df['skip']==False].apply(lambda x: background_ratio(x['position'],x['strand'],\
        ko_datasets[x['strand']],w=w,g=g,winsorization_level=winsorization_level,mode=mode), axis=1)

    ref_df.loc[ref_df['skip']==False, 'ref_bg_depth'] = ref_df.loc[ref_df['skip']==False].apply(lambda x: background_ratio(x['position'],x['strand'],\
        ref_datasets[x['strand']],w=w,g=g,winsorization_level=winsorization_level, mode=mode, return_depth=True), axis=1)
    ko_df.loc[ko_df['skip']==False, 'ko_bg_depth'] = ko_df.loc[ko_df['skip']==False].apply(lambda x: background_ratio(x['position'],x['strand'],\
        ko_datasets[x['strand']],w=w,g=g,winsorization_level=winsorization_level, mode=mode, return_depth=True), axis=1)

    ref_df = ref_df[['position','strand','ref_bg_ratio', 'ref_bg_depth']]
    ko_df = ko_df[['position','strand','ko_bg_ratio', 'ko_bg_depth']]

    out_df = ref_df.merge(ko_df, on=['position','strand'], how='outer')
    out_df['ko_id'] = ko_id
    # pdb.set_trace()
    return out_df

def get_background_ratios(ref_em, ko_ems, ko_ids, countcol = 'counts', w=50, g=3, winsorization_level=0.0, skip_positions=None, mode = 'up'):
    print('Calculating all background ratios...')
    out_df = [calc_background_ratio(ref_em, ko_ems[i], ko_ids[i], w=w, g=g, winsorization_level = winsorization_level, countcol=countcol, mode=mode, skip_positions=skip_positions) for i in range(len(ko_ids))]
    print('Concatenating background tables...')
    return pd.concat(out_df)

#region functions to determine if there is a paired 3' peak
def get_threeprime_peak_ratio(fiveprime_em, threeprime_em, w, g, skip_positions = None, winsor_level = 0, predicted_offset=1,countcol = 'counts', outcol = 'minus1_3end_bg_ratio', mode='up', keep_skip=True):
    print('Calculating threeprime flags...')

    fiveprime_df = fiveprime_em.readcounts.reset_index()
    threeprime_df = threeprime_em.readcounts.reset_index()

    # Skip positions is defined as the positions with signal that exceed a certain threshold
    fiveprime_df.loc[fiveprime_df['strand']=='+','skip'] = skip_positions['+']
    fiveprime_df.loc[fiveprime_df['strand']=='-','skip'] = skip_positions['-']
    threeprime_df.loc[threeprime_df['strand']=='+','skip'] = skip_positions['+']
    threeprime_df.loc[threeprime_df['strand']=='-','skip'] = skip_positions['-']

    threeprime_datasets = {'+': threeprime_df.loc[threeprime_df['strand']=='+',countcol].to_list(),'-': threeprime_df.loc[threeprime_df['strand']=='-',countcol].to_list()}

    fiveprime_df.loc[fiveprime_df['skip']==False,outcol] = fiveprime_df.loc[fiveprime_df['skip']==False].apply(
        lambda x: background_ratio(x['position'], x['strand'], threeprime_datasets[x['strand']],w=w,g=g, winsorization_level=winsor_level, mode=mode, p_offset=predicted_offset),axis=1)

    if keep_skip:
        fiveprime_df = fiveprime_df[['position','strand','skip',outcol]]
    else:
        fiveprime_df = fiveprime_df[['position','strand',outcol]]

    return fiveprime_df
#endregion


def calc_zscore(ref_em, ko_em, ko_id, z_depth_thresh=0, skip_positions=None, countcol = 'counts', w=50, g=3, winsorization_level=0.0, use_winsor=False):
    print('Calculating z-scores...')

    ref_df = ref_em.readcounts.reset_index()
    ko_df = ko_em.readcounts.reset_index()

    # Skip positions is defined as the positions with signal that exceed a certain threshold
    ref_df.loc[ref_df['strand']=='+','skip'] = skip_positions['+']
    ref_df.loc[ref_df['strand']=='-','skip'] = skip_positions['-']
    ko_df.loc[ko_df['strand']=='+','skip'] = skip_positions['+']
    ko_df.loc[ko_df['strand']=='-','skip'] = skip_positions['-']

    ref_datasets = {'+': ref_df.loc[ref_df['strand']=='+',countcol].to_list(),'-': ref_df.loc[ref_df['strand']=='-',countcol].to_list()}
    ko_datasets = {'+': ko_df.loc[ko_df['strand']=='+',countcol].to_list(),'-': ko_df.loc[ko_df['strand']=='-',countcol].to_list()}

    ref_df.loc[ref_df['skip']==False, 'ref_zscore'] = ref_df.loc[ref_df['skip']==False].apply(lambda x: zscore(x['position']-1,w,\
        ref_datasets[x['strand']],g, z_depth_thresh, winsor_thresh=winsorization_level, local_winsorize=use_winsor), axis=1)
    ko_df.loc[ko_df['skip']==False, 'ko_zscore'] = ko_df.loc[ko_df['skip']==False].apply(lambda x: zscore(x['position']-1,w,\
        ko_datasets[x['strand']],g, z_depth_thresh, winsor_thresh=winsorization_level, local_winsorize=use_winsor), axis=1)

    ref_df = ref_df[['position','strand','ref_zscore']]
    ko_df = ko_df[['position','strand','ko_zscore']]

    out_df = ref_df.merge(ko_df, on=['position','strand'], how='outer')
    out_df['ko_id'] = ko_id
    return out_df    

def get_zscores(ref_em, ko_ems, ko_ids, countcol = 'counts', w=50, g=3, z_depth_thresh=0, winsorization_level=0.0, use_winsor=False, skip_positions=None):
    print('Calculating all zscores...')
    out_df = [calc_zscore(ref_em, ko_ems[i], ko_ids[i], w=w, g=g, z_depth_thresh=z_depth_thresh, winsorization_level = winsorization_level, countcol=countcol, skip_positions=skip_positions) for i in range(len(ko_ids))]
    print('Concatenating zscore tables...')
    return pd.concat(out_df)
# endregion

# region Functions for grouping together adjacent positions (which can be grouped using groupby)
def get_groups(df, metric, threshold, id_to_use, max_gap = 1):
    print('Getting groups...')
    # To simplify, digitize data
    df = df[df['ko_id']==id_to_use].copy(deep=True)
    df.loc[df[metric]>=threshold,'passed'] = 1
    df.loc[df[metric]<threshold,'passed'] = 0
    df.loc[np.isnan(df[metric]),'passed'] = 0  # Many may be nan, due to skipped positions in metric calculation

    pass_bools = {'+': df.loc[df['strand']=='+', 'passed'].to_list(),'-': df.loc[df['strand']=='-','passed'].to_list()}

    # Add column specifying distance to next value passing threshold
    next_distances_plus = np.zeros(len(pass_bools['+']))
    for i in range(len(pass_bools['+'])):
        for j in range(1,max_gap+2):
            if i+j>len(pass_bools['+'])-1:
                next_distances_plus[i] = np.nan
                break  
            downstream_bool = pass_bools['+'][i+j]
            if downstream_bool == 1:
                next_distances_plus[i] = j
                break
            if j>max_gap+1:  # Don't want to waste time iterating past max gap
                next_distances_plus[i] = np.nan
                break
            next_distances_plus[i] = np.nan  # Catch positions at the end of pass_bools 
        
    next_distances_minus = np.zeros(len(pass_bools['-']))
    for i in range(len(pass_bools['-'])):
        for j in range(1,max_gap+2):
            if i+j>len(pass_bools['-'])-1:
                next_distances_minus[i] = np.nan
                break  
            downstream_bool = pass_bools['-'][i+j]
            if downstream_bool == 1:
                next_distances_minus[i] = j
                break
            if j>max_gap+1:  # Don't want to waste time iterating far past max gap
                next_distances_minus[i] = np.nan
                break
            next_distances_minus[i] = np.nan  # Catch positions at the end of pass_bools 

    df.loc[df['strand']=='+', 'bool_dist'] = next_distances_plus
    df.loc[df['strand']=='-', 'bool_dist'] = next_distances_minus

    plus_groups = []
    in_group = False
    group_n = 0
    for i,d in enumerate(next_distances_plus):
        if pass_bools['+'][i] == 0:
            if np.isnan(d) or d>max_gap:
                if not in_group:
                    plus_groups.append(group_n)
                    group_n+=1
                else:
                    group_n+=1
                    plus_groups.append(group_n)
                    group_n+=1
                    in_group = False
            elif (d <= max_gap) and ~np.isnan(d) and in_group:
                plus_groups.append(group_n)
            elif (d <= max_gap) and ~np.isnan(d) and ~in_group:
                plus_groups.append(group_n)
                group_n+=1
        elif pass_bools['+'][i] == 1:
            plus_groups.append(group_n)
            in_group = True

    minus_groups = []
    in_group = False
    group_n+=1
    for i,d in enumerate(next_distances_minus):
        if pass_bools['-'][i] == 0:
            if np.isnan(d) or d>max_gap:  # Positions outside of group
                if not in_group:
                    minus_groups.append(group_n)
                    group_n+=1
                else:  # Will hit this just downstream of group
                    group_n+=1 # Added to differentiate position just downstream of group
                    minus_groups.append(group_n)
                    group_n+=1
                    in_group = False
            elif (d <= max_gap) and ~np.isnan(d) and in_group: # Positions in gaps
                minus_groups.append(group_n)
            elif (d <= max_gap) and ~np.isnan(d) and ~in_group: # Position
                minus_groups.append(group_n)
                group_n+=1
        elif pass_bools['-'][i] == 1:
            minus_groups.append(group_n)
            in_group = True

    df.loc[df['strand']=='+','group'] = plus_groups
    df.loc[df['strand']=='-','group'] = minus_groups
    df = df[['position','strand','group']]

    return df
# endregion

# region Final function to call all the above
def build_endseq_table(ref_5end_em, ko_5end_ems, ref_rend5_em, ko_rend5_ems, ref_rend3_em, ref_3end_em, ko_ids, genome_seq, group_thresh, gap=2, half_window=50, 
                    window_3end=50, gap_3end=2, winsor=0.05, data_column = 'cds_normalized', pseudocount_val = 0.01,
                     predicted_offset=1, skip_thresh=10, winsor_zscore = False, z_depth_thresh=0, group_metric='ref_bg_ratio'):

    # Parameters for depth/step calling
    skip_pos = get_skip_positions(ref_5end_em, skip_thresh)  # Could use these positions to filter l

    rat = get_all_ratios(ref_5end_em, ko_5end_ems, ko_ids, countcol = data_column, ref_suffix = '_ref_5end', ko_suffix = '_ko_5end')

    rend = get_depth_step_rend_filters(ref_rend5_em, ko_rend5_ems, ko_ids, half_window, gap, skip_positions=skip_pos, winsor_level = winsor, countcol = data_column, mode='down')
    three = get_threeprime_peak_ratio(ref_5end_em, ref_3end_em, w=window_3end, g=gap_3end, skip_positions = skip_pos, winsor_level = winsor, countcol = data_column, predicted_offset=predicted_offset, mode='up')
    three_minus3 = get_threeprime_peak_ratio(ref_5end_em, ref_3end_em, w=window_3end, g=gap_3end, skip_positions = skip_pos, winsor_level = winsor, \
        countcol = data_column, predicted_offset=3, keep_skip=False, outcol='minus3_3end_bg_ratio', mode='up')
    zsc = get_zscores(ref_5end_em, ko_5end_ems, ko_ids, countcol=data_column, z_depth_thresh=z_depth_thresh, w=half_window, g = gap, winsorization_level=winsor, use_winsor=winsor_zscore, skip_positions=skip_pos)
    bg = get_background_ratios(ref_5end_em, ko_5end_ems, ko_ids, countcol=data_column, w=half_window, g=gap, winsorization_level=winsor, skip_positions=skip_pos,mode='down')
    group = get_groups(df=bg, metric=group_metric, threshold=group_thresh, id_to_use=ko_ids[0],max_gap=1) # This would need to use a different df if a different metric is used (e.g. z-score)
    kmer = get_site_kmers(ref_5end_em, 20, genome_seq, skip_positions=skip_pos)
    dz = get_deadzone(ref_rend3_em, w=50, zthresh=8, d_thresh=3, gap=2, av_thresh=0.25, em_type='3end', deadzone_dist=40, local_winsorize=False, skip_thresh=3, winsor_thresh = 0.05, save_z = False, outcol='deadzone_5endseq')

    a = rend.merge(three, on=['position','strand'], how='outer')
    b = a.merge(zsc, on=['position','strand','ko_id'], how='outer')
    c = b.merge(bg, on=['position','strand','ko_id'], how='outer')
    d = c.merge(rat, on=['position','strand','ko_id'], how='outer')
    e = d.merge(group, on=['position','strand'], how='outer')
    f = e.merge(three_minus3, on=['position','strand'], how='outer')
    g = f.merge(dz, on=['position','strand'], how='outer')
    outdf = g.merge(kmer, on=['position', 'strand'], how='outer')

    outdf['ref_pseudocounted'] = False
    outdf['ko_pseudocounted'] = False

    outdf.loc[outdf[data_column+'_ref_5end']==0.0,'ref_pseudocounted'] = True
    outdf.loc[outdf[data_column+'_ref_5end']==0.0,data_column+'_ref_5end'] = pseudocount_val

    outdf.loc[outdf[data_column+'_ko_5end']==0.0,'ko_pseudocounted'] = True
    outdf.loc[outdf[data_column+'_ko_5end']==0.0,data_column+'_ko_5end'] = pseudocount_val

    outdf.loc[outdf['5end_ratio']!=outdf['5end_ratio'],'5end_ratio'] = outdf[data_column+'_ref_5end']/pseudocount_val

    return outdf
#endregion

# region Write sequences to files

def write_sequences_from_lists(list_dict, w, genome, outF, fasta=False):
    with open(outF, 'w') as f:
        for p in list_dict['+']:
            seq = get_sequences_from_position(p, '+', w, genome)
            if fasta:
                f.write('>%s_%s\n' % (seq, str(p)))
            f.write(seq+'\n')
        for p in list_dict['-']:
            seq = get_sequences_from_position(p, '-', w, genome)
            if fasta:
                f.write('>%s_%s\n' % (seq, str(p)))
            f.write(seq+'\n')

def write_sequences_from_df(peak_df, w, genome, outF, fasta=False):
    f_peaks = peak_df.loc[peak_df['strand']=='+','position'].to_list()
    r_peaks = peak_df.loc[peak_df['strand']=='-','position'].to_list()
    with open(outF, 'w') as f:
        for p in f_peaks:
            seq = get_sequences_from_position(p, '+', w, genome)
            if fasta:
                f.write('>%s_%s\n' % (seq, str(p)))
            f.write(seq+'\n')
        for p in r_peaks:
            seq = get_sequences_from_position(p, '-', w, genome)
            if fasta:
                f.write('>%s_%s\n' % (seq, str(p)))
            f.write(seq+'\n')

def write_kplogo_meme_from_df(out_base, peak_df, half_length, genome_seq):
    write_sequences_from_df(peak_df, half_length, genome_seq, out_base+'_sequences.fa', fasta=True)
    write_sequences_from_df(peak_df, half_length, genome_seq, out_base+'_sequences.txt', fasta=False)

def write_kplogo_meme_from_lists(out_base, peak_lists, half_length, genome_seq):
    write_sequences_from_lists(peak_lists, half_length, genome_seq, out_base+'_sequences.fa', fasta=True)
    write_sequences_from_lists(peak_lists, half_length, genome_seq, out_base+'_sequences.txt', fasta=False)

#endregion

def strand_specific_groupby(df, count_col = 'counts'):

    plus_grouped = df.loc[df['strand']=='+']
    plus_grouped['ko_group_sum'] = plus_grouped.groupby(['group'])['%s_ko_5end' % count_col].transform('sum')
    plus_grouped['ref_group_sum'] = plus_grouped.groupby(['group'])['%s_ref_5end' % count_col].transform('sum')
    plus_grouped['pos_in_group'] = plus_grouped.groupby(['group'])['position'].transform('count') # Count number of positions within the group

    plus_grouped['ko_pseudocounted'] = plus_grouped['pos_in_group'] - \
        plus_grouped.groupby(['group'])['ko_pseudocounted'].transform('sum')   # if all positions in ko are pseudocounted, this will be zero
    plus_grouped['ko_pseudocounted'] = plus_grouped['ko_pseudocounted']==0  # When all positions in group pseudocounted, flag it as True

    # assign to max position in (-) strand, min in (+) strand
    min_position = plus_grouped[["group","position"]].groupby("group").min().reset_index()
    plus_grouped = pd.merge(plus_grouped, min_position, on = ['position', 'group'], how="inner")

    minus_grouped = df.loc[df['strand']=='-']
    minus_grouped['ko_group_sum'] = minus_grouped.groupby(['group'])['%s_ko_5end' % count_col].transform('sum')
    minus_grouped['ref_group_sum'] = minus_grouped.groupby(['group'])['%s_ref_5end' % count_col].transform('sum')
    minus_grouped['pos_in_group'] = minus_grouped.groupby(['group'])['position'].transform('count') # Count number of positions within the group

    minus_grouped['ko_pseudocounted'] = minus_grouped['pos_in_group'] - \
        minus_grouped.groupby(['group'])['ko_pseudocounted'].transform('sum')   # if all positions in ko are pseudocounted, this will be zero
    minus_grouped['ko_pseudocounted'] = minus_grouped['ko_pseudocounted']==0  # When all positions in group pseudocounted, flag it as True

    # assign to max position in (-) strand, min in (+) strand
    max_position = minus_grouped[["group","position"]].groupby("group").max().reset_index()
    minus_grouped = pd.merge(minus_grouped, max_position, on = ['position', 'group'], how="inner")

    df_grouped = pd.concat([plus_grouped,minus_grouped],axis=0)
    return df_grouped

def strand_specific_groupby_maxthresh(df, count_col = 'counts'):
    '''
    Group positions together which were flagged as "grouped" due to adjacency and exceeding peak-to-bg ratio thresh
    Calculate the max value of peak-to-background thresholds within these windows to be used in later thresholding
    This is done to avoid skipping positions in the middle of peaks when calculating the sum over that peak
    '''

    plus_grouped = df.loc[df['strand']=='+']
    plus_grouped['ko_group_sum'] = plus_grouped.groupby(['group'])['%s_ko_5end' % count_col].transform('sum')
    plus_grouped['ref_group_sum'] = plus_grouped.groupby(['group'])['%s_ref_5end' % count_col].transform('sum')
    plus_grouped['max_ref_bg_ratio'] = plus_grouped.groupby(['group'])['ref_bg_ratio'].transform('max')
    plus_grouped['max_ko_bg_ratio'] = plus_grouped.groupby(['group'])['ko_bg_ratio'].transform('max')
    plus_grouped['max_minus1_3end_bg_ratio'] = plus_grouped.groupby(['group'])['minus1_3end_bg_ratio'].transform('max')
    plus_grouped['avg_ref_depth'] = plus_grouped.groupby(['group'])['ref_depth'].transform(np.mean)
    plus_grouped['avg_ko_depth'] = plus_grouped.groupby(['group'])['ko_depth'].transform(np.mean)

    min_position = minus_grouped[["group","position"]].groupby("group").min().reset_index()
    plus_grouped = pd.merge(plus_grouped, min_position, on = ['position', 'group'], how="inner")
    
    minus_grouped = df.loc[df['strand']=='-']
    minus_grouped['ko_group_sum'] = minus_grouped.groupby(['group'])['%s_ko_5end' % count_col].transform('sum')
    minus_grouped['ref_group_sum'] = minus_grouped.groupby(['group'])['%s_ref_5end' % count_col].transform('sum')
    minus_grouped['max_ref_bg_ratio'] = minus_grouped.groupby(['group'])['ref_bg_ratio'].transform('max')
    minus_grouped['max_ko_bg_ratio'] = minus_grouped.groupby(['group'])['ko_bg_ratio'].transform('max')
    minus_grouped['max_minus1_3end_bg_ratio'] = minus_grouped.groupby(['group'])['minus1_3end_bg_ratio'].transform('max')
    minus_grouped['avg_ref_depth'] = minus_grouped.groupby(['group'])['ref_depth'].transform(np.mean)
    minus_grouped['avg_ko_depth'] = minus_grouped.groupby(['group'])['ko_depth'].transform(np.mean)

    max_position = minus_grouped[["group","position"]].groupby("group").max().reset_index()
    minus_grouped = pd.merge(minus_grouped, max_position, on = ['position', 'group'], how="inner")
    
    df_grouped = pd.concat([plus_grouped,minus_grouped],axis=0)
    return df_grouped

def get_cleavage_positions(df, ref_bg_ratio_thresh=7.5, ref_depth_thresh=0.5, minus1_3end_bg_ratio_thresh=7.5, ndoA_sites=False):
    if not ndoA_sites:
        return df[(df['ref_bg_ratio']>ref_bg_ratio_thresh) & (df['ref_depth']>ref_depth_thresh) & (df['minus1_3end_bg_ratio']>minus1_3end_bg_ratio_thresh) & (df['deadzone_5endseq']==False)]
    else:
        return df[(df['ref_bg_ratio']>ref_bg_ratio_thresh) & (df['ref_depth']>ref_depth_thresh) & (df['minus3_3end_bg_ratio']>minus1_3end_bg_ratio_thresh) & (df['deadzone_5endseq']==False)]

def get_cleavage_positions_maxthresh(df, ref_bg_ratio_thresh=7.5, minus1_3end_bg_ratio_thresh=7.5, ref_depth_thresh=0.5, ko_depth_thresh=0.5, ndoA_sites=False):
    if not ndoA_sites:
        return df[(df['max_ref_bg_ratio']>ref_bg_ratio_thresh) & (df['avg_ref_depth']>ref_depth_thresh) \
                  & (df['max_minus1_3end_bg_ratio']>minus1_3end_bg_ratio_thresh) & (df['avg_ko_depth']>ko_depth_thresh)]
    else:
        return df[(df['max_ref_bg_ratio']>ref_bg_ratio_thresh) & (df['avg_ref_depth']>ref_depth_thresh) \
                  & (df['max_minus3_3end_bg_ratio']>minus1_3end_bg_ratio_thresh) & (df['avg_ko_depth']>ko_depth_thresh)]







