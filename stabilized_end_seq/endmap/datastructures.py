import pandas as pd
import numpy as np
import pdb
import matplotlib.pyplot as plt
from plotting.utils import cdf
import pickle
from rend_utilities.wig_utils import find_nearby_values, calculate_step, calculate_local_depth, check_nearby_peak, gene_from_position, sum_region
from rend_utilities.peak_parsing_utils import is_in_coding, is_near_coding
from Bio.Seq import Seq
from annotation_utils.expression_quantification import find_none_genes
from peak_calling.zscore import calculate_all_z, group_nearby_peaks
from scipy.stats.mstats import winsorize

class EndMap(object):
    '''
    Stores read counts, imported from a pair of .wig files and associated tables of
    peaks called from these arrays of read counts.

    Attributes
    ----------
    type : str
        Either '3end' or '5end'. This may not currently be used.
    readcounts : pandas.DataFrame
        Dataframe containing read counts, typically imported directly from .wig files
    peaks : tuple of lists
        Lists of forward and reverse strand peaks, called with method call_peaks_zscore
    ID : str
        Used to store an ID associated with a dataset
    ee_table : pandas.DataFrame
        If calculating and importing end enrichment table using JBL scripts, store here
    peaktable : pandas.DataFrame
        Contains peak locations and other associate information about them (e.g. steps)
        Calculated by method call_peaks_zscore

    Methods
    -------
    read_wigs(data_dir, forward_wiggle, reverse_wiggle, count_colname = 'counts', \
        chrom='NC_000964.3', chrom_length=4215606)
        Reads a pair of .wigs into self.readcounts attribute
    calculate_RPM()
        Calculates reads per million column ('RPM') and adds to self.readcounts
    normalize_to_cds(mochi_ref, colname = 'cds_normalized')
        Calculates a column in self.readcounts where counts at each position are
        normalized to the total counts within all coding regions. Calculates this 
        by reading a mochiview gene location set (mochi_ref) and summing over these regions
    plot_zscore_distribution(outf, w, gap=0, av_thresh=0, local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05)
        Calculates the z scores at each position and generates a CDF of the distribution of all z-scores
    z_sweep(outf, w, gap, av_thresh,local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05)
        Calculates and orders the z scores at each position and then plots log2 z-scores for all positions
        Should ideally see a clear elbow in this distribution for defining z-score cutoff
    call_peaks_zscore(w, zthresh, d_thresh, gap, av_thresh, local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05, save_z = False)
        Calls positions that exceed a particular Z-score threshold (zthresh)
        Saves results in self.peaks and self.peaktable
    write_peak_positions(outf)
        Writes .csv of peak positions
    calculate_peak_steps(window, gap, data_colname = 'counts', col_name = 'step',winsorization_level=0)
        Calculates step in read density associated with each peak in self.peaktable and saves as new column
    calculate_peak_depth(window, gap, data_colname = 'counts', col_name = 'depth',winsorization_level=0)
        Calculates local read depth (winsorized mean) around peaks in self.peaktable and saves as new column

    
    '''

    def __init__(self):
        self.type = None # '3end' or '5end'
        self.readcounts = None # Dataframe with read counts read from .wig files
        self.peaks = None
        self.ID = None

        self.ee_table = None # Table of "end enrichment" file generated as output of Jean's peak calling scripts
        self.peaktable = None # dataframe containing peaks associated w/ dataset

    def read_wigs(self, data_dir, forward_wiggle, reverse_wiggle, count_colname = 'counts', chrom='NC_000964.3', chrom_length=4215606):

        df_f = pd.read_csv(data_dir+forward_wiggle,delimiter='\t',header=None,names=['position',count_colname],dtype={'position': int, count_colname:float})
        df_r = pd.read_csv(data_dir+reverse_wiggle,delimiter='\t',header=None,names=['position',count_colname],dtype={'position': int, count_colname:float})

        df_f['strand'] = '+'
        df_r['strand'] = '-'

        df_f.set_index('position',inplace=True)
        df_r.set_index('position',inplace=True)

        zero_counts = pd.DataFrame({'position':np.array(range(1,chrom_length+1)),count_colname:np.zeros(chrom_length)})
        zero_counts.set_index('position',inplace=True)

        full_df_f = zero_counts.merge(df_f,on='position',how='left')
        full_df_f.fillna(0.0,inplace=True)
        full_df_f[count_colname] = full_df_f['%s_x' % (count_colname)] + full_df_f['%s_y' % (count_colname)]
        full_df_f=full_df_f[[count_colname]]
        full_df_f['strand'] = '+'

        full_df_r = zero_counts.merge(df_r,on='position',how='left')
        full_df_r.fillna(0.0,inplace=True)
        full_df_r[count_colname] = full_df_r['%s_x' % (count_colname)] + full_df_r['%s_y' % (count_colname)]
        full_df_r=full_df_r[[count_colname]]
        full_df_r['strand'] = '-'

        self.readcounts = pd.concat([full_df_f,full_df_r])

    def calculate_RPM(self):
        assert self.readcounts is not None
        self.readcounts['RPM'] = (self.readcounts['counts']/sum(self.readcounts['counts']))*1000000

    def normalize_to_cds(self, mochi_ref, colname = 'cds_normalized'):
        assert self.readcounts is not None
        mochi = pd.read_csv(mochi_ref,delimiter='\t')
        # Filter to real genes only
        mochi['has_none_gene'] = mochi['FEATURE_NAME'].apply(find_none_genes) 
        cds = mochi[(mochi['ALIASES']=='gene') & (mochi['has_none_gene']==False)]

        f_counts = self.readcounts.loc[self.readcounts['strand']=='+','counts']
        r_counts = self.readcounts.loc[self.readcounts['strand']=='-','counts']

        cds['reads'] = cds.apply(lambda x: sum_region(f_counts, r_counts, x['START'], x['END'], x['STRAND']),axis=1)

        self.readcounts[colname] = (self.readcounts['counts']/sum(cds['reads']))*1000000

    def normalize_to_cds_winsor(self, mochi_ref, winsorization, colname = 'cds_normalized'):
        assert self.readcounts is not None
        mochi = pd.read_csv(mochi_ref,delimiter='\t')
        # Filter to real genes only
        mochi['has_none_gene'] = mochi['FEATURE_NAME'].apply(find_none_genes) 
        cds = mochi[(mochi['ALIASES']=='gene') & (mochi['has_none_gene']==False)]

        f_counts = self.readcounts.loc[self.readcounts['strand']=='+','counts']
        r_counts = self.readcounts.loc[self.readcounts['strand']=='-','counts']

        cds['reads'] = cds.apply(lambda x: sum_region(f_counts, r_counts, x['START'], x['END'], x['STRAND']),axis=1)

        self.readcounts[colname] = (self.readcounts['counts']/sum(winsorize(cds['reads'].to_numpy(), winsorization)))*1000000
        
    def plot_zscore_distribution(self, outf, w, gap=0, av_thresh=0, local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05):
        data_f = self.readcounts.loc[self.readcounts['strand']=='+','counts'].to_list()
        z_array_f = calculate_all_z(w,data_f,gap=gap,av_thresh=av_thresh,local_winsorize=local_winsorize,winsor_thresh=winsor_thresh,skip_thresh=skip_thresh)
        data_r = self.readcounts.loc[self.readcounts['strand']=='-','counts'].to_list()
        z_array_r = calculate_all_z(w,data_r,gap=gap,av_thresh=av_thresh,local_winsorize=local_winsorize,winsor_thresh=winsor_thresh,skip_thresh=skip_thresh)

        z_array_f_noNA = list(z_array_f[np.where((z_array_f == z_array_f) & (~np.isinf(z_array_f)))])
        z_array_r_noNA = list(z_array_r[np.where((z_array_r == z_array_r) & (~np.isinf(z_array_r)))])

        with open('z_array_f_noNA', 'wb') as f:
            pickle.dump(z_array_f_noNA,f)
        with open('z_array_r_noNA', 'wb') as f:
            pickle.dump(z_array_r_noNA,f)

        cdf(z_array_f_noNA+z_array_r_noNA, flipy=True)
        # plt.semilogy()
        plt.xlim([0,1])
        plt.ylim([-20,20])
        plt.savefig(outf,format='png',dpi=250)
        plt.close()

    def z_sweep(self, outf, w, gap, av_thresh,local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05):
        data_f = self.readcounts.loc[self.readcounts['strand']=='+','counts'].to_list()
        z_array_f = calculate_all_z(w,data_f,gap,av_thresh,local_winsorize=local_winsorize,winsor_thresh=winsor_thresh,skip_thresh=skip_thresh)
        data_r = self.readcounts.loc[self.readcounts['strand']=='-','counts'].to_list()
        z_array_r = calculate_all_z(w,data_r,gap,av_thresh,local_winsorize=local_winsorize,winsor_thresh=winsor_thresh,skip_thresh=skip_thresh)

        # Drop Z scores which have values which are <1, np.nan, or np.inf
        z_array_f_noNA = list(z_array_f[np.where((z_array_f == z_array_f) & (z_array_f>1) & (~np.isinf(z_array_f)))])
        z_array_r_noNA = list(z_array_r[np.where((z_array_r == z_array_r) & (z_array_r>1) & (~np.isinf(z_array_r)))])

        # Consider both strands together
        all_z_list = z_array_f_noNA+z_array_r_noNA

        x = np.arange(len(all_z_list))
        y = sorted(np.log2(all_z_list))

        plt.plot(x,y)
        plt.axhline(y=3)
        plt.ylim([0,10])
        plt.savefig(outf,format='png',dpi=250)
        plt.close()

    def call_peaks_zscore(self, w, zthresh, d_thresh, gap, av_thresh, local_winsorize=False, skip_thresh=0, winsor_thresh = 0.05, save_z = False):
        print('Calling peaks...')
        data_f = self.readcounts.loc[self.readcounts['strand']=='+','counts'].to_list()
        z_array_f = calculate_all_z(w,data_f,gap, av_thresh, local_winsorize=local_winsorize,winsor_thresh=winsor_thresh,skip_thresh=skip_thresh)
        f_peaks = np.where(z_array_f>zthresh)[0]  # [0] required as np.where returns a tuple, but this is a 1D array
        f_peak_z = z_array_f[np.where(z_array_f>zthresh)]

        data_r = self.readcounts.loc[self.readcounts['strand']=='-','counts'].to_list()
        z_array_r = calculate_all_z(w,data_r,gap, av_thresh, local_winsorize=local_winsorize,winsor_thresh=winsor_thresh,skip_thresh=skip_thresh)
        r_peaks = np.where(z_array_r>zthresh)[0]  # [0] required as np.where returns a tuple, but this is a 1D array
        r_peak_z = z_array_r[np.where(z_array_r>zthresh)]

        if save_z:
            self.readcounts.loc[self.readcounts['strand']=='+','z_score'] = z_array_f
            self.readcounts.loc[self.readcounts['strand']=='-','z_score'] = z_array_r

        # Group together positions which are within d_thresh of one another. Take position with highest z-score (or first position if identical)
        print('Grouping peaks...')
        f_peaks, f_peak_z = group_nearby_peaks(f_peaks, f_peak_z, d_thresh)
        r_peaks, r_peak_z = group_nearby_peaks(r_peaks, r_peak_z, d_thresh)
        # pdb.set_trace()

        # Switch over to 1-based indexing, equivalent to what the end-enrichment files from JBL scripts output
        f_peaks = f_peaks + 1
        r_peaks = r_peaks + 1

        self.peaks = (f_peaks,r_peaks)  # Lists of peaks stored as tuple in self.peaks
        # pdb.set_trace()
        pt_dict = {'position':list(f_peaks)+list(r_peaks),
                   'strand':['+']*len(f_peaks)+['-']*len(r_peaks),
                   'zscore':list(f_peak_z)+list(r_peak_z)}

        self.peaktable = pd.DataFrame(pt_dict)  # For future operations DataFrame of peaks stored

    def write_peak_positions(self,outf):
        self.peaktable.to_csv(outf,sep='\t')

    def calculate_peak_steps(self, window, gap, data_colname = 'counts', col_name = 'step',winsorization_level=0):
        assert self.peaktable is not None and self.readcounts is not None
        df = self.peaktable
        dsets = {'+': self.readcounts.loc[self.readcounts['strand']=='+',data_colname].to_list(),
         '-': self.readcounts.loc[self.readcounts['strand']=='-',data_colname].to_list()}

        df[col_name] = df.apply(lambda x: calculate_step(x['position'], x['strand'], dsets[x['strand']],window,gap,winsorization_level=winsorization_level),axis=1)

    def calculate_peak_depth(self, window, gap, data_colname = 'counts', col_name = 'depth',winsorization_level=0):
        assert self.peaktable is not None and self.readcounts is not None
        df = self.peaktable
        dsets = {'+': self.readcounts.loc[self.readcounts['strand']=='+',data_colname].to_list(),
         '-': self.readcounts.loc[self.readcounts['strand']=='-',data_colname].to_list()}

        df[col_name] = df.apply(lambda x: calculate_local_depth(x['position'], x['strand'], dsets[x['strand']],window,gap,winsorization_level=winsorization_level),axis=1)
