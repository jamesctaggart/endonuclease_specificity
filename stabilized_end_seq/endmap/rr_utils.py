from scipy.optimize import curve_fit
from scipy.stats import gamma, skewnorm
from endmap.calc_endseq_columns import get_cleavage_positions, strand_specific_groupby, get_cleavage_positions_maxthresh, strand_specific_groupby_maxthresh
from rend_utilities.wig_utils import find_nearby_values
import numpy as np
import pandas as pd
import pdb
import matplotlib.pyplot as plt
import pickle
from brokenaxes import brokenaxes

def get_cleavage_sites_RR(es_df, RR_thresh, ref_counts_thresh, ko_id, ref_depth_thresh, \
    ko_depth_thresh, ref_bg_ratio_thresh, minus1_3end_bg_ratio_thresh, ppnorm, data_column, ndoA_sites=False):

    summary_df = es_df.loc[(es_df['ko_id']==ko_id) & (es_df['ko_depth']>ko_depth_thresh)]
    putative_sites = get_cleavage_positions(summary_df, ref_bg_ratio_thresh=ref_bg_ratio_thresh, \
        minus1_3end_bg_ratio_thresh=minus1_3end_bg_ratio_thresh,ref_depth_thresh=ref_depth_thresh,ndoA_sites=ndoA_sites)
    gps = strand_specific_groupby(putative_sites, count_col=data_column)
    if not ppnorm:
        gps['ko_normalized_peak_height'] = gps['ko_group_sum']/gps['ko_depth']
        gps['ref_normalized_peak_height'] = gps['ref_group_sum']/gps['ref_depth']
    else:     # If ppnorm, normalize to total reads mapping to peak pairs rather than total reads mapping to coding regions
        gps['ko_normalized_peak_height'] = (gps['ko_group_sum']/sum(gps['ko_group_sum']))/gps['ko_depth']
        gps['ref_normalized_peak_height'] = (gps['ref_group_sum']/sum(gps['ref_group_sum']))/gps['ref_depth']

    gps['DB_RR'] = gps['ref_normalized_peak_height']/gps['ko_normalized_peak_height']

    
    if data_column == 'counts': # This threshold of <1 read may not work if not using raw counts, so skip if not
        gps.loc[gps['ko_group_sum']<1,'pseudocounted'] = True
        gps.loc[gps['ko_group_sum']>=1,'pseudocounted'] = False

    gps = gps.loc[gps['ref_group_sum']>=ref_counts_thresh]
    sites = gps.loc[(gps['DB_RR']>RR_thresh)]

    # gps_noKO_filter provides all grouped sites without thresholding on a particular KO dataset
    # If you take union of peaks in gps_noKO_filter from all datasets, will give maximal number of sites we might be able to assign
    no_ko_filter_sites = get_cleavage_positions(es_df.loc[(es_df['ko_id']==ko_id)], ref_bg_ratio_thresh=ref_bg_ratio_thresh, \
        minus1_3end_bg_ratio_thresh=minus1_3end_bg_ratio_thresh,ref_depth_thresh=ref_depth_thresh,ndoA_sites=ndoA_sites)
    gps_noKO_filter = strand_specific_groupby(no_ko_filter_sites, count_col=data_column)
    gps_noKO_filter = gps_noKO_filter.loc[gps_noKO_filter['ref_group_sum']>=ref_counts_thresh]

    return sites, gps, putative_sites, summary_df, gps_noKO_filter

def get_cleavage_sites_ko_bg_ratio(es_df, ko_bg_ratio_thresh, ref_counts_thresh, ko_id, ref_depth_thresh, \
    ko_depth_thresh, ref_bg_ratio_thresh, minus1_3end_bg_ratio_thresh, data_column='counts',ndoA_sites=False):

    summary_df = es_df.loc[(es_df['ko_id']==ko_id) & (es_df['ko_depth']>ko_depth_thresh)]
    putative_sites = get_cleavage_positions(summary_df, ref_bg_ratio_thresh=ref_bg_ratio_thresh, \
        minus1_3end_bg_ratio_thresh=minus1_3end_bg_ratio_thresh,ref_depth_thresh=ref_depth_thresh,ndoA_sites=ndoA_sites)
    gps = strand_specific_groupby(putative_sites, count_col=data_column)
    sites = gps.loc[(gps['ko_bg_ratio']<ko_bg_ratio_thresh) & (gps[data_column+'_ref_5end']>=ref_counts_thresh)]
    return sites, gps, putative_sites,summary_df

def get_RR_background_distribution(es_df, ref_counts_thresh, ko_id, ref_depth_thresh, \
    ko_depth_thresh, data_column='counts'):

    df = es_df.loc[(es_df['ko_id']==ko_id) & (es_df['ko_depth']>ko_depth_thresh) & (es_df['ref_depth']>ref_depth_thresh)]
    df = df.loc[(df[data_column+'_ref_5end']>ref_counts_thresh) & (df[data_column+'_ko_5end']>=1)]
    df['ko_normalized_peak_height'] = df[data_column+'_ko_5end']/df['ko_depth']
    df['ref_normalized_peak_height'] = df[data_column+'_ref_5end']/df['ref_depth']
    df['DB_RR'] = df['ref_normalized_peak_height']/df['ko_normalized_peak_height']
    return np.log10(np.array(df['DB_RR']))

def gaussian(x, mu, sigma, A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def gaussian_sum(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return gaussian(x,mu1,sigma1,A1)+gaussian(x,mu2,sigma2,A2)

def fit_two_gaussian(data, initial_params, ref_dist = None, outbase='gauss_fit'):
    print('Fitting %s...' % outbase)
    assert len(initial_params) == 6
    y,x,_=plt.hist(data,100,alpha=.3,label='data')
    x=(x[1:]+x[:-1])/2 # make sure that len(x)==len(y)

    params,cov=curve_fit(gaussian_sum,x,y,initial_params)
    sigma=np.sqrt(np.diag(cov))
    plt.close()
    return params, sigma

def plot_2gauss_fit(data, params, outbase='gauss_fit', vline=None, xlim=[-2,5],ylim=None, get_bins=False, pseudo_count=None, pseudoX=None):

    _,x,_=plt.hist(data,100,alpha=1,label='data',color='k',zorder=2)
    x=(x[1:]+x[:-1])/2 # enforce len(x)==len(y)

    fit_x_vals = np.linspace(x.min(), x.max(), 600)

    plt.plot(fit_x_vals, gaussian_sum(fit_x_vals, *params), color='red', lw=2, label='Sum',zorder=3)
    plt.plot(fit_x_vals, gaussian(fit_x_vals, *params[:3]), color='red', lw=1, ls="--", label='Gaussian 1',zorder=3)
    plt.plot(fit_x_vals, gaussian(fit_x_vals, *params[3:]), color='red', lw=1, ls="--", label='Gaussian 2',zorder=3)

    if vline is not None:
        plt.axvline(x=vline,c='k',linestyle='--',linewidth=1.5)

    if pseudo_count is not None and pseudoX is not None:
        width = x[1]-x[0]
        extended_bins = width*np.arange(-100,300)
        plt.bar(pseudoX,pseudo_count,width=width,label='pseudo_vals',color='gray') 
    plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.savefig('figures/%s.pdf' % outbase, transparent=True)
    plt.close()
    
    if get_bins:
        return x


def flag_pseudocounted(gps, normalized=False, grouped = True):

    if grouped:
        if not normalized:
            gps.loc[gps['ko_group_sum']<1,'pseudocounted'] = True
            gps.loc[gps['ko_group_sum']>=1,'pseudocounted'] = False
        else:
            gps.loc[gps['ko_group_sum']<0.01,'pseudocounted'] = True
            gps.loc[gps['ko_group_sum']>=0.01,'pseudocounted'] = False   
    else:
        if not normalized:
            gps.loc[gps['counts_ko_5end']<1,'pseudocounted'] = True
            gps.loc[gps['counts_ko_5end']>=1,'pseudocounted'] = False
        else:
            gps.loc[gps['counts_ko_5end']<0.01,'pseudocounted'] = True
            gps.loc[gps['counts_ko_5end']>=0.01,'pseudocounted'] = False  

    return gps

def get_non_pseudo_RR(gps, normalized=False, grouped = True):
    gps = flag_pseudocounted(gps, normalized=normalized, grouped=grouped)
      
    return np.log10(gps.loc[gps['pseudocounted']==False,'DB_RR'].to_numpy())

def bimodal_gaussian_params_to_dict(params):
    return {'mu1':params[0],'sigma1':params[1], 'a1':params[2],'mu2':params[3],'sigma2':params[4],'a2':params[5]}


def count_sites_multidf(gps1, gps2, dist_thresh, save_sites=False, site_file=''):
    plus_pos_1 = gps1.loc[gps1['strand']=='+'].position.to_list()
    plus_pos_2 = gps2.loc[gps2['strand']=='+'].position.to_list()
    minus_pos_1 = gps1.loc[gps1['strand']=='-'].position.to_list()
    minus_pos_2 = gps2.loc[gps2['strand']=='-'].position.to_list()

    plus_nearby = find_nearby_values(plus_pos_1,plus_pos_2,dist_thresh=dist_thresh, symmetric=True)
    minus_nearby = find_nearby_values(minus_pos_1,minus_pos_2,dist_thresh=dist_thresh, symmetric=True)

    # Count the number of sites which are shared with distance 0<x≤dist_thresh
    nearby_count = 0
    for nearby_list in plus_nearby, minus_nearby:
        for pair in nearby_list:
            if pair[0] == pair[1]:
                continue
            else:
                nearby_count+=1

    # Get the union of the two dataframes
    merged = gps1.merge(gps2,how='outer',on=['position','strand'], indicator=True)

    if save_sites:
        with open(site_file, 'wb') as f:
            pickle.dump(merged, f)

    return len(merged)-nearby_count

    


