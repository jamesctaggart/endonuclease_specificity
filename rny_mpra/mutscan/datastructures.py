import pickle
import pandas as pd
import numpy as np
from mutscan.importutils import map_umi_to_vbc, quantify_barcodes, extract_variant_sequences
from mutscan.sequtils import get_mutations, count_mutations, make_all_mutations, \
    make_all_double_mutations, find_mutations, combine_mutation_strings, check_mutation_string, \
    offset_mutations
from mutscan.heatmap_utils import nest_dm_dict
from mutscan.bc_correction import calculate_bccorr
from mutscan.structutils import fold_all
from Bio import SeqIO
from collections import OrderedDict
from mutscan.plotting import snp_boxplot, get_fig_axes, snp_boxplot_topbars
import seaborn as sns
import matplotlib.pyplot as plt
import pdb
import os
import copy
from collections import defaultdict

class MutScanPool(object):
    def __init__(self, unmutated_seq, zero_position):
        self.unmutated_seq = unmutated_seq  # Length of this should be length of gDNA BC-to-seq read
        self.read_length = len(unmutated_seq)
        self.zero_position = zero_position # 5' end position

        self.data = None # Dataframe containing all info about variants / counts
        self.bc_to_seq = None  # dictionary {BC: variant_sequence}
        self.wt_median = None
        self.mutated_region = None
        self.df_pickle_path = ''
        self.RNAquant_pickle_path = ''
        self.DNAquant_pickle_path = ''
        self.bcmap_pickle_path = ''
        self.raw_path = ''
        self.threshold = None  # Can use this to record threshold used in threshold_data()
        self.parsing_positions = {'umi_start':None, 'umi_end':None, 'bc_start':None, 'bc_end':None,
                              'constant_region_start':None, 'constant_region_end':None}
        self.constant_region_seq = None
        self.depth_threshold = None
        self.identifier = None  # Short string for ID when writing file names, etc
        self.double_mutant_df = None # Store all-background (neutral only) all by all mutations

    def read_from_raw(self, RNA_bc_fastq, gDNA_bc_fastq, gDNA_seq_fastq1, gDNA_seq_fastq2):

        RNA_bc_read = SeqIO.parse(RNA_bc_fastq, "fastq")
        gDNA_bc_read = SeqIO.parse(gDNA_bc_fastq, "fastq")
        gDNA_seq_read1 = SeqIO.parse(gDNA_seq_fastq1, "fastq")
        gDNA_seq_read2 = SeqIO.parse(gDNA_seq_fastq2, "fastq")


        RNA_counts = quantify_barcodes(RNA_bc_read, umi_start = self.parsing_positions['umi_start'],
                                       umi_end=self.parsing_positions['umi_end'],
                                       vbc_start = self.parsing_positions['bc_start'],
                                       vbc_end=self.parsing_positions['bc_end'],
                                       constant_region_start=self.parsing_positions['constant_region_start'],
                                       constant_region_end=self.parsing_positions['constant_region_end'],
                                       constant_region_seq = self.constant_region_seq)

        gDNA_counts = quantify_barcodes(gDNA_bc_read, umi_start = self.parsing_positions['umi_start'],
                                       umi_end=self.parsing_positions['umi_end'],
                                       vbc_start = self.parsing_positions['bc_start'],
                                       vbc_end=self.parsing_positions['bc_end'],
                                       constant_region_start=self.parsing_positions['constant_region_start'],
                                       constant_region_end=self.parsing_positions['constant_region_end'],
                                       constant_region_seq = self.constant_region_seq)

        self.bc_to_seq = extract_variant_sequences(gDNA_seq_read1, gDNA_seq_read2, umi_start = self.parsing_positions['umi_start'],
                                       umi_end=self.parsing_positions['umi_end'],
                                       vbc_start = self.parsing_positions['bc_start'],
                                       vbc_end=self.parsing_positions['bc_end'],
                                       constant_region_start=self.parsing_positions['constant_region_start'],
                                       constant_region_end=self.parsing_positions['constant_region_end'],
                                       constant_region_seq = self.constant_region_seq)

        self.data = initialize_variant_dataframe(RNA_counts, gDNA_counts, self.bc_to_seq, self.unmutated_seq)

    def read_from_pickle(self, bcmap_fh, rnaquant_fh, dnaquant_fh, df_fh = '', precomputed_df = False):
        with open(self.bcmap_pickle_path+bcmap_fh,'rb') as f:
            # print(f)
            self.bc_to_seq = pickle.load(f)

        if precomputed_df:
            with open(self.df_pickle_path+df_fh,'rb') as f:
                self.data = pickle.load(f)
        else:
            with open(self.RNAquant_pickle_path+rnaquant_fh, 'rb') as f:
                RNA_counts = pickle.load(f)
            with open(self.DNAquant_pickle_path+dnaquant_fh, 'rb') as f:
                gDNA_counts = pickle.load(f)

            self.data = initialize_variant_dataframe(RNA_counts, gDNA_counts, self.bc_to_seq, self.unmutated_seq)

    def write_pickle(self, fh_base):
        with open(self.pickle_path+fh_base+'_dataframe','w') as f:
            pickle.dump(self.data, f)

        with open(self.pickle_path+fh_base+'_BCtoseq','w') as f:
            pickle.dump(self.bc_to_seq, f)

    def report_bc_bias(self, out_fh, minimum_variants=10, max_mutations=2):
        print('Writing normalized barcode biases...')
        outf = open(out_fh,'w')

        df = self.data.loc[self.data.n_mutations<=max_mutations,:]

        for seq in df.sequence.unique():
            seq_rows = df.loc[df.sequence == seq,['barcode','RNA:gDNA']]
            if len(seq_rows)>=minimum_variants:
                print(seq,len(seq_rows))
                seq_median = np.median(seq_rows['RNA:gDNA'])
                for bc in seq_rows.barcode:
                    normalized_bc_ratio = float(seq_rows.loc[seq_rows.barcode == bc,'RNA:gDNA']/seq_median)
                    outf.write(str(bc)+'\t'+str(normalized_bc_ratio)+'\n')
        print('Finished writing barcode biases.')
        outf.close()

    def merge_pool(self, pool_to_merge_in):
        '''
        Concatenates the dataframes contained in MutScanPool.data 

        :param pool_to_merge_in: A MutScanPool object
        :return:
        '''
        self.data = pd.concat([self.data,pool_to_merge_in.data])

    def normalize_data(self):
        df = self.data
        self.wt_median = np.median(df[df['sequence'] == self.unmutated_seq]['RNA:gDNA'])
        df['norm_RNA:gDNA'] = df['RNA:gDNA'] / self.wt_median

    def threshold_data(self, threshold, column):
        df = self.data
        self.threshold = str(threshold)
        self.data = df[df[column]>threshold]

    def calculate_bc_normalization(self):
        bc_corr = calculate_bccorr(self.data, self.unmutated_seq)
        self.data['bc_correction'] = self.data['barcode'].apply(lambda x: bc_corr[x])
        self.data['corrected_ratio'] = self.data.apply(lambda x: x['norm_RNA:gDNA']/bc_corr[x.barcode],axis=1)

    def call_mutations(self):
        df = self.data
        df['mutations'] = df['sequence'].apply(get_mutations, consensus_seq=self.unmutated_seq)
        df['n_mutations'] = df['sequence'].apply(count_mutations, consensus_seq=self.unmutated_seq)

    def get_snps(self, bc_threshold = None,include_wt=False, column = 'corrected_ratio', start_idx=None, end_idx=None):
        snp_dict = {}

        # Call mutations if not yet called
        try:
            self.data['mutations']
        except(KeyError):
            self.call_mutations()

        # Return dictionary of all single-nucleotide mutations
        mut_list = make_all_mutations(self.unmutated_seq, start_idx=start_idx, end_idx=end_idx)
        for m in mut_list:
            snp_dict[m] = self.data.loc[self.data['mutations']==m, column].to_list()

        if bc_threshold is not None:
            for m in mut_list:
                if len(snp_dict[m])<bc_threshold:
                    snp_dict.pop(m)

        if include_wt: # Changed so it doesn't require .unmutated_seq to match.
            snp_dict['wt'] = self.data.loc[self.data['mutations']=='', column].to_list()
        # if include_wt:
        #     snp_dict['wt'] = self.data.loc[self.data['sequence']==self.unmutated_seq, column].to_list()

        return snp_dict

    def get_snps_allbg(self, bc_threshold = None, include_wt=False, start_idx=None, end_idx=None,
                       only_neutral=True,neutral_threshold=0.9, column = 'corrected_ratio'):
        snp_dict = {}

        # Call mutations if not yet called
        try:
            self.data['mutations']
        except(KeyError):
            self.call_mutations()

        if only_neutral:
            flagged_snps = self.get_nonneutral_snps(neutral_threshold, column=column)

        # Return dictionary of all single-nucleotide mutations
        mut_list = make_all_mutations(self.unmutated_seq, start_idx=start_idx, end_idx=end_idx)
        for m in mut_list:
            if only_neutral:
                mut_idx = self.data.apply(lambda x: find_mutations(x['mutations'],m, mutations_to_avoid=flagged_snps),axis=1)
            else:
                mut_idx = self.data.apply(lambda x: find_mutations(x['mutations'],m),axis=1)

            snp_dict[m] = self.data.loc[mut_idx, column].to_list()


        if bc_threshold is not None:
            for m in mut_list:
                if len(snp_dict[m])<bc_threshold:
                    snp_dict.pop(m)

        if include_wt:
            snp_dict['wt'] = self.data.loc[self.data['sequence']==self.unmutated_seq, column].to_list()

        return snp_dict

    def get_neutral_snps(self, bc_threshold = None, threshold = 0.9, column='corrected_ratio'):
        snps = self.get_snps(bc_threshold = bc_threshold, include_wt=True, column=column)
        snps = {k:np.median(v) for k,v in snps.items()}
        neutral_snps = {k:v for k,v in snps.items() if abs(np.log2(v/snps['wt']))<abs(np.log2(threshold))}
        neutral_snps.pop('wt')  # WT will always be flagged as neutral, so remove
        return neutral_snps.keys()

    def get_nonneutral_snps(self, bc_threshold = None, threshold = 0.9, column='corrected_ratio'):
        snps = self.get_snps(bc_threshold = bc_threshold, include_wt=True, column=column)
        snps = {k:np.median(v) for k,v in snps.items()}  # Collapse all barcodes to median
        flagged_snps = {k:v for k,v in snps.items() if abs(np.log2(v/snps['wt']))>abs(np.log2(threshold))}
        return flagged_snps.keys()

    def flag_nonneutral(self, neutral_threshold=0.9, column='corrected_ratio'):
        '''

        :param neutral_threshold: Threshold (in log2) above/below which a snp will get flagged
        :param column: Column to evaluate "neutrality" over - typically some RNA/DNA ratio for this MRPA
        :return: None. Adds column to self.data.
        '''
        # Add a column to flag variants which contain mutations which are not neutral in a single mutant context

        self.data['contains_nonneutral_mutation'] = False

        flagged_snps = self.get_nonneutral_snps(neutral_threshold, column=column)

        for snp in flagged_snps:
            idx = self.data.apply(lambda x: find_mutations(x['mutations'], snp), axis=1)
            self.data.loc[idx,'contains_nonneutral_mutation'] = True

    def make_snp_boxplot(self, axes=None, bc_threshold = None, allbg = False, start_idx=None, end_idx=None,
                         only_neutral=True,neutral_threshold=0.9,logy=False, column ='corrected_ratio', show_wt=True,
                         hlines=[],figsize=[30,10],topbar=True,record_zeros=None, figname='figures/snp_boxplot.png', **kwargs):
        if allbg:
            snps = self.get_snps_allbg(bc_threshold=bc_threshold, include_wt=True, column=column, start_idx=start_idx, end_idx=end_idx,
                                       only_neutral=only_neutral, neutral_threshold=neutral_threshold)
        else:
            snps = self.get_snps(bc_threshold=bc_threshold, include_wt=True, column=column, start_idx=start_idx, end_idx=end_idx)

        if record_zeros is not None:
            record_zeros.write(self.identifier+'\n')

        snp_boxplot(snps,logy=logy, show_wt=show_wt,hlines=hlines,figsize=figsize,figname=figname,topbar=topbar,record_zeros=record_zeros,**kwargs)


    def make_conditional_snp_boxplot(self, conditioning_mutations, axes=None, bc_threshold=None, allbg=False, start_idx=None, end_idx=None,
                         only_neutral=True, neutral_threshold=0.9, logy=False, column='corrected_ratio',
                         show_wt=True,hlines=[], **kwargs):
        if allbg:
            snps = self.get_conditioned_snps_allbg(conditioning_mutations, bc_threshold=bc_threshold, include_wt=True, column=column, start_idx=start_idx,
                                       end_idx=end_idx,
                                       only_neutral=only_neutral, neutral_threshold=neutral_threshold)
        else:
            snps = self.get_conditioned_snps(conditioning_mutations ,bc_threshold=bc_threshold, include_wt=True, column=column, start_idx=start_idx, end_idx=end_idx)

        fig, ax = get_fig_axes(axes)

        snp_boxplot(snps,logy=logy, axes=ax,show_wt=show_wt,hlines=hlines, **kwargs)

    def get_conditioned_snps(self,conditioning_mutations, bc_threshold = None, include_wt=False, include_cond_nosnp=True, column = 'corrected_ratio', start_idx=None, end_idx=None):
        '''

        :param conditioning_seq:  Sequence used as the starting point to make mutations from. Replaces self.unmutated
        :param bc_threshold:
        :param include_wt:
        :param include_cond_nosnp:
        :param column:
        :param start_idx:
        :param end_idx:
        :return:
        '''
        # Get the values associated w/ each point mutation conditioned on a non-wildtype sequence
        # Can use this to look at the impacts of second snps on a first snp

        snp_dict = {}

        # Call mutations if not yet called
        try:
            self.data['mutations']
        except(KeyError):
            self.call_mutations()

        # Return dictionary of all single-nucleotide mutations
        mut_list = make_all_mutations(self.unmutated_seq, start_idx=start_idx, end_idx=end_idx)
        for m in mut_list:
            if find_mutations(conditioning_mutations,m):  # Skip mutations which are in the condition
                continue
            try:# Add in conditioning mutations to each
                m_combined = combine_mutation_strings(m,conditioning_mutations)
            except(AssertionError):
                continue
            snp_dict[m] = self.data.loc[self.data['mutations']==m_combined, column].to_list()

        if bc_threshold is not None:
            for m in mut_list:
                if len(snp_dict[m])<bc_threshold:
                    snp_dict.pop(m)

        if include_wt:
            snp_dict['wt'] = self.data.loc[self.data['sequence']==self.unmutated_seq, column].to_list()
        if include_cond_nosnp:
            snp_dict['cond_mut_%s' % conditioning_mutations] = self.data.loc[self.data['mutations']==conditioning_mutations, column].to_list()

        return snp_dict

    def get_conditioned_snps_allbg(self, conditioning_mutations, bc_threshold = None, include_wt=False,include_cond_nosnp=True,
                                   start_idx=None, end_idx=None, only_neutral=True,neutral_threshold=0.9, column = 'corrected_ratio'):

        snp_dict = {}

        # Call mutations if not yet called
        try:
            self.data['mutations']
        except(KeyError):
            self.call_mutations()

        if only_neutral:
            flagged_snps = self.get_nonneutral_snps(neutral_threshold, column=column)

        # Return dictionary of all single-nucleotide mutations
        mut_list = make_all_mutations(self.unmutated_seq, start_idx=start_idx, end_idx=end_idx)
        for m in mut_list:
            print(m)

            try:# Add in conditioning mutations to each
                m = combine_mutation_strings(m,conditioning_mutations)
            except(AssertionError): # Skip cases where two mutations at same position (which violates an assert in fn)
                continue

            if only_neutral:
                mut_idx = self.data.apply(lambda x: find_mutations(x['mutations'],m, mutations_to_avoid=flagged_snps),axis=1)
            else:
                mut_idx = self.data.apply(lambda x: find_mutations(x['mutations'],m),axis=1)

            snp_dict[m] = self.data.loc[mut_idx, column].to_list()


        if bc_threshold is not None:
            for m in mut_list:
                if len(snp_dict[m])<bc_threshold:
                    snp_dict.pop(m)

        if include_wt:
            snp_dict['wt'] = self.data.loc[self.data['sequence']==self.unmutated_seq, column].to_list()
        if include_cond_nosnp:
            snp_dict['cond_mut_%s' % conditioning_mutations] = self.data.loc[self.data['mutations']==conditioning_mutations, column].to_list()

        return snp_dict

    def make_condensed_mutdict(self, filter_nonneutral=False, nmut_thresh = 2, column='corrected_ratio'):
        '''
        This function is used to generate a plotting-ready dictionary of values associated with each mutation string
        Parameters allow pre-filtering of this dictionary.

        Used my functions such as get_double_mutants (used in making 2-mutant heatmaps).

        :param filter_nonneutral: Set true if you want to exclude variants which have >nmut_thresh mutations, 1+ of which is not neutral
        :param nmut_thresh: Number of mutations above which you filter out non-neutral mutations
        :param column: Column to return in value of dictionary
        :return: dict {mutation_string : list of values from df with that mutation, ready for plotting}
        '''

        df = self.data
        out = dict()

        if filter_nonneutral:# For double mutant analysis, nmut_thresh = 2. For single mutant, =1.
            subdf = df.loc[(df['contains_nonneutral_mutation']==False) | (df['n_mutations']<=nmut_thresh)]
        else:
            subdf = df

        for m in subdf['mutations'].unique():
            out[m] = list(subdf.loc[subdf['mutations']==m, column])

        return out

    def get_double_mutants(self, start_idx=None, end_idx=None, bc_threshold = None, include_wt=False,
                           allbg=False, only_neutral=True, column='corrected_ratio'):
        '''
        Extracts lists of all mutant impacts from dataframe using make_condensed_mutdict() and then finds
        the value in specified column for all variants w/ each possible double mutant

        Used to extract values to plot in make_double_mutant_plots()

        :param start_idx: Start position to extract dat
        :param end_idx:
        :param bc_threshold: Mutation strings which have fewer than bc_threshold variants will be dropped if not None
        :param include_wt: If your plot needs a 'wt' set to true
        :param allbg: If true, include variants with additional mutations as long as they have primary 2
        :param only_neutral: If true and allbg==True, only include variants lacking "non-neutral" background mutations -> neutral previously calculated by self.flag_nonneutral()
        :param column: Column containing value over which to plot.
        :return:
        '''

        double_mutants = make_all_double_mutations(self.unmutated_seq, start_idx=start_idx, end_idx=end_idx)

        print('Double mutant plots: condensing data...')
        if allbg:
            if only_neutral:  # If only_neutral, remove non-neutral-containing mutations when you condense data to dict
                condensed_data = self.make_condensed_mutdict(filter_nonneutral = True, column=column)
            else:
                condensed_data = self.make_condensed_mutdict(filter_nonneutral = False, column=column)
        else:
            condensed_data = self.make_condensed_mutdict(filter_nonneutral = True, column=column)

        dm_dict = {}
        print('Double mutant plots: collecting ratios...')
        for dm in double_mutants:
            dm_dict[dm] = []
            print(dm)
            if allbg:
                for key in condensed_data:
                    if find_mutations(key, dm):
                        dm_dict[dm] += condensed_data[key]
            else:
                try:
                    dm_dict[dm] = condensed_data[dm]
                except(KeyError):
                    print('Missing in double mutants: '+dm)


            if bc_threshold is not None:
                if len(dm_dict[dm])<bc_threshold:
                    dm_dict.pop(dm)

        if include_wt:
            dm_dict['wt'] = self.data.loc[self.data['mutations']=='', column].to_list()

        return dm_dict

    def get_single_mutant_products(self, start_idx=None, end_idx=None, column='corrected_ratio'):
        single_mutants = make_all_mutations(self.unmutated_seq, start_idx=start_idx, end_idx=end_idx)

        singleproduct_df = {}
        for m1 in single_mutants:
            singleproduct_df[m1] = {}
            for m2 in single_mutants:
                singleproduct_df[m1][m2] = \
                    np.log2(np.median(
                        self.data[self.data['mutations'] == m1][column]) * \
                            np.median(
                                self.data[self.data['mutations'] == m2][column]))
        return singleproduct_df

    def calculate_double_mutant_dataframe(self, start_idx=None, end_idx=None, column='corrected_ratio',
                                 allbg=True, only_neutral=True,log2_average=True):


        # Extract values associated w/ each double mutant
        doublemut_dict = self.get_double_mutants(start_idx=start_idx, end_idx=end_idx, column = column, allbg=allbg, only_neutral=only_neutral)
        doublemut_dict = nest_dm_dict(doublemut_dict,log2_average=log2_average)  # LOG2 MEDIAN CALCULATED HERE
        doublemut_df = pd.DataFrame(doublemut_dict)

        self.double_mutant_df = doublemut_df

    def save_double_mutant_plot(self, out_name, x_start=None,x_end=None,y_start=None,y_end=None,
                                vmin=None,vmax=None,figx=38,figy=30,flipy=False,flipx=False):
        plt.figure(figsize=(figx,figy))
        df = self.double_mutant_df
        if flipy:
            df = df.loc[make_all_mutations(self.unmutated_seq,x_start,x_end)][::-1][make_all_mutations(self.unmutated_seq,y_start,y_end)]
        else:
            df = df.loc[make_all_mutations(self.unmutated_seq,x_start,x_end)][make_all_mutations(self.unmutated_seq,y_start,y_end)]

        if flipx:
            df = df.iloc[:, ::-1]

        mask = df.isnull()
        sns.heatmap(df, mask=mask, vmin=vmin, vmax=vmax,
                    cmap='RdYlBu')  # For color see: https://www.r-graph-gallery.com/38-rcolorbrewers-palettes
        plt.savefig(out_name, format='pdf', dpi=300)
        plt.close()

    def calculate_RNA_secondary_structure(self, start_idx, end_idx):
        # First, write all the sequences to a file
        # pdb.set_trace()
        with open('temp_rna_seqs','w') as f:
            for s in self.data['sequence'].to_list():
                f.write(s[start_idx:end_idx]+'\n')

        # Then feed the file containing all the sequences into RNAfold, extracting the MFE and associated dot-bracket secondary structure.
        dG_array, fold_array = fold_all('temp_rna_seqs')
        self.data['deltaG'] = dG_array
        self.data['fold'] = fold_array

        os.remove('temp_rna_seqs')

    def get_collapsed_variants(self):
        # Return a dataframe for which all variants with the same sequence are collapsed to one
        # Creates new columns for the number of barcodes with that sequence, as well as mean/median ∆G and normalized RNA:gDNA
        gb = self.data.groupby(['sequence'])
        counts = gb.size().to_frame(name='barcode_count')
        return (counts\
        .join(gb.agg({'norm_RNA:gDNA': 'mean'}).rename(columns={'norm_RNA:gDNA': 'mean_norm_RNA:gDNA'}))\
        .join(gb.agg({'norm_RNA:gDNA': 'median'}).rename(columns={'norm_RNA:gDNA': 'median_norm_RNA:gDNA'}))\
        .join(gb.agg({'deltaG': 'mean'}).rename(columns={'deltaG': 'mean_deltaG'}))\
        .join(gb.agg({'deltaG': 'median'}).rename(columns={'deltaG': 'median_deltaG'}))\
        .reset_index())

    def get_mut_filtered_data(self, muts_of_interest, only_neutral=True, 
                              clean_bg = True, neutral_threshold=0.9, column='corrected_ratio'):
        '''
        Use this to isolate subsets of the dataframe which contain predefined sets of mutations.
        This can be used for coarse-grained analysis of variants which of mutations of a certain class (e.g. in stem or UA-rich downstream element)
        '''

        flagged_snps = self.get_nonneutral_snps(neutral_threshold, column=column)

        self.data['has_%s_neutralbg' % muts_of_interest] = self.data.apply(lambda x: check_mutation_string(x['mutations'], muts_of_interest, flagged_snps), axis=1)
        
        # pdb.set_trace()

        return self.data.loc[self.data['has_%s_neutralbg' % muts_of_interest]==True]
    
    def get_filtered_on_mut_list(self, mutations_list, only_neutral=True, 
                              clean_bg = True, neutral_threshold=0.9, column='corrected_ratio'):
        '''
        Use this to isolate subsets of the dataframe which contain predefined sets of mutations.
        This can be used for coarse-grained analysis of variants which of mutations of a certain class (e.g. in stem or UA-rich downstream element)
        '''

        df_list = []
        for m in mutations_list:
            print(m)
            df_list.append(self.get_mut_filtered_data(m, neutral_threshold=neutral_threshold, column=column))
        # pdb.set_trace()
        return pd.concat(df_list).drop_duplicates()

    def threshold_sweep(self):
        x = self.data.sort_values('gDNA_count')

        binsize=40
        bins = [(binsize*i,binsize*(i+1)) for i in range(0,20)]

        data_lists = []
        stdevs = []
        for bin in bins:

            bin_data = x[(x['sequence']==self.unmutated_seq) & (x['gDNA_count']>=bin[0]) & 
                         (x['gDNA_count']<bin[1])]['RNA:gDNA'].to_numpy()
            stdevs.append(np.std(bin_data))
            data_lists.append(bin_data)

            print(bin, np.std(bin_data), len(bin_data))

        boxprops = dict(linestyle='-', linewidth=1, color='black')
        medianprops = dict(linestyle='-', linewidth=3, color='k')
        whiskerprops = dict(linestyle='--', color='black')

        plt.boxplot([np.log2(x) for x in data_lists], whis=[5, 95],
        showfliers=False,medianprops=medianprops,boxprops=boxprops,whiskerprops=whiskerprops, patch_artist=True,zorder=3)
        plt.savefig('figures/parameter_tuning/%s_binned_data.pdf' % self.identifier, format='pdf')
        plt.close()

        plt.plot([b[0] for b in bins],stdevs)
        plt.xlabel('Minimum reads in bin (binsize: %s)' % str(binsize))
        plt.ylabel('Standard Deviation')
        plt.ylim([0,4])
        plt.savefig('figures/parameter_tuning/%s_threshold_vs_stdev.pdf' % self.identifier, format='pdf')
        plt.close()

        plt.plot([b[0] for b in bins],[np.log10(len(d)) for d in data_lists])
        plt.xlabel('Minimum reads in bin (binsize: %s)' % str(binsize))
        plt.ylabel('Log10 number of datapoints')
        # plt.semilogy()
        plt.savefig('figures/parameter_tuning/%s_threshold_vs_count_datapoints.pdf' % self.identifier, format='pdf')
        plt.close()

def initialize_variant_dataframe(RNA_quantification, gDNA_quantification, gDNA_variants, unmutated_seq):

    '''
    :RNA_quantificiation: Dictionary of form {VBC: number of UMIs for variant}
    :gDNA_quantificiation: Dictionary of form {VBC: number of UMIs for variant}
    :gDNA_variants: Dictionary of form {VBC: variant sequence}
    :unmutated_variable_region: Reference sequence of unmutated variable region
    '''

    # Pull out the list of RNA/DNA ratios for each variant.
    # Each number in the list associated w/ a variant is a different barcoded sequence
    # Only focus on variants which exceed depth threshold to avoid noise in ratios
    variant_ratios = list()
    variant_gDNA_counts = list()
    variant_RNA_counts = list()
    variant_barcodes = list()
    variants = list()
    missing_in_full_gDNA = set()

    for vbc in gDNA_quantification:
        ratio = RNA_quantification[vbc]/float(gDNA_quantification[vbc])
        if gDNA_quantification[vbc]:
            try:
                x,y,z = gDNA_variants[vbc], RNA_quantification[vbc], gDNA_quantification[vbc] # This is to trigger KeyError
                variant_ratios.append(ratio)
                variants.append(str(gDNA_variants[vbc]))
                variant_RNA_counts.append(RNA_quantification[vbc])
                variant_gDNA_counts.append(gDNA_quantification[vbc])
                variant_barcodes.append(vbc)
            except(KeyError):
                missing_in_full_gDNA.add(vbc)

    outdf = pd.DataFrame({'sequence':variants, 'barcode':variant_barcodes,'gDNA_count':variant_gDNA_counts,
                          'RNA_count':variant_RNA_counts, 'RNA:gDNA':variant_ratios})
    # outdf.set_index('sequence',inplace=True)

    print('missing in gDNA ',len(missing_in_full_gDNA))

    return outdf


def get_mutation_values(msp,mutation_list, neutralbg=True):
    if neutralbg:
        nonneutral_snps = msp.get_nonneutral_snps(0.9, column='norm_RNA:gDNA')
    else:
        nonneutral_snps = []

    mutation_values = {}
    included_mutations = {}

    for m in mutation_list:
        idx = msp.data.apply(lambda x: find_mutations(x['mutations'], m, mutations_to_avoid=nonneutral_snps), axis=1)
        mutation_values[m] = msp.data.loc[idx, 'norm_RNA:gDNA'].values
        included_mutations[m] = msp.data.loc[idx, 'mutations'].tolist()

    return mutation_values, included_mutations

def merge_pools(msplist, idlist, new_unmutated_sequence, unmutated_seq_offsets):
    '''
    Use MutScanPool.merge_pool() to merge together multiple MSP objects.
    This is useful when plotting upstream and downstream scanning together!

    unmutated_seq_offsets is a list of int that is the offset that needs to be applied
    to mutation numbers. This enables SNP plots, etc.
    '''
    count_thresh = msplist[0].count_threshold
    for i,msp in enumerate(msplist):
        msp.data['expt_id'] = idlist[i]
        msp.data['seq_index_offset'] = unmutated_seq_offsets[i]
        assert count_thresh==msp.count_threshold
    
    # Initialize the output. It will lack attribute definitions, need to populate these.
    merged_msp = MutScanPool(new_unmutated_sequence, 0)
    merged_msp.identifier = msplist[0].identifier
    merged_msp.unmutated_seq = new_unmutated_sequence
    merged_msp.data = msplist[0].data

    # Merge together the data table for all included MutScanPools
    for msp in msplist[1:]:
        merged_msp.identifier = merged_msp.identifier+'---'+msp.identifier
        merged_msp.merge_pool(msp)

    # Use the provided list of offsets to correct the mutation string for each variant
    merged_msp.data['mutations'] = merged_msp.data.apply(
        lambda x: offset_mutations(x['mutations'],x['seq_index_offset']),axis=1)
    merged_msp.count_threshold = count_thresh
    
    return merged_msp
