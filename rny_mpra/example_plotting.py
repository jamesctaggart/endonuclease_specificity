from mutscan.datastructures import merge_pools
from mutscan.sequtils import get_mutations, count_mutations
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib
matplotlib.style.use('jct_light_spline')


records = open('not_plotted_counts.txt', 'w')  # Record the number of variants which have 0 RNA counts or are outside bounds of the plots.

offsets = [0,len('GCTTACTTTAAAAAGCcacgcaacacggttctcgtc')]  # Showing cggR sequence here, change for other sites
extended_cggR_seq = 'GCTTACTTTAAAAAGCcacgcaacacggttctcgtcacagacgaaggagccgcaaagaagttattaagggatgaataatccctcaatataaatatctctcacttatttaaa'.upper()
# endregion

records.write('Number of variants with a zero ratio in each dataset:'+'\n')

# region Generate SNP plot for the full cggR region (FL detection)
msp_list = [msp_upstream_mutations, msp_downstream_mutations] # These are can be initialized as in example_loading_msp.py
msp_merged = merge_pools(msp_list, [msp.identifier for msp in msp_list], extended_cggR_seq, offsets)  # Merge together all comparable experiments targeting the same cleavage site 

figsize=[30,10]; ylim=[-2,4.5]; figname='%s_thresh:%s_cleanbg_logy.pdf'%(msp_merged.identifier,str(msp_merged.count_threshold))
msp_merged.make_snp_boxplot(figname=figname, figsize=figsize, ylim=ylim, column='norm_RNA:gDNA',allbg=False, logy=True, start_idx=30, end_idx=95, show_wt=False,
                            vert_stripe=False, topbar=True, hlines=[-2,-1,0,1,2,3,4], record_zeros = records)


# Generate double mutant plots for all variants
precomputed = True
for msp in [msp_merged]:  # Iterate through all experiments, if multiple being considered
    msp.df_pickle_path = 'pickles/'
    if not precomputed:
        msp.calculate_double_mutant_dataframe(start_idx=30, end_idx=95,column='norm_RNA:gDNA')
        with open(msp.df_pickle_path+'_'+msp.identifier+'_'+str(msp.count_threshold)+'_doublemutdf','wb') as f:
            pickle.dump(msp.double_mutant_df, f)
    else:
        with open(msp.df_pickle_path+'_'+msp.identifier+'_'+str(msp.count_threshold)+'_doublemutdf','rb') as f:
            msp.double_mutant_df = pickle.load(f)

x_start=30; x_end=95; y_start=30; y_end=95; figx=9; figy=7
msp_merged.save_double_mutant_plot('%s_neutralbg_%sthresh_doublemutant_%s-%sby%s-%s.pdf' %
                            (msp_merged.identifier,str(msp_merged.count_threshold),str(x_start),str(x_end),str(y_start),str(y_end)),
                            x_start=x_start,x_end=x_end,y_start=y_start,y_end=y_end,vmin=-3,vmax=3,figx=figx,figy=figy)


# Generate secondary structure plots for all variants
min_barcodes = 3
markers = ['o','v',',']
cthresh=msp_upstream_mutations.count_threshold # Can define threshold for number of min counts per variant to be included in plot 

colors = ['black']
fig, ax = plt.subplots()
y_limits_FL = [-2,4.5]
out_of_bounds = 0
for i, msp in enumerate([msp_downstream_mutations]):

    df = msp.get_collapsed_variants()

    unmutated_dG = df.loc[df['sequence']==msp.unmutated_seq]['median_deltaG'].values[0]

    collapsed = msp.get_collapsed_variants()
    collapsed['mutations'] = collapsed['sequence'].apply(get_mutations, consensus_seq=msp.unmutated_seq)
    collapsed['n_mutations'] = collapsed['sequence'].apply(count_mutations, consensus_seq=msp.unmutated_seq)

    collapsed_filtered = collapsed.loc[(collapsed['barcode_count']>=min_barcodes) & (collapsed['n_mutations']>0)]

    mut1 = collapsed_filtered.loc[collapsed_filtered['n_mutations']==1]
    mut2 = collapsed_filtered.loc[collapsed_filtered['n_mutations']==2]
    mut3plus = collapsed_filtered.loc[collapsed_filtered['n_mutations']>=3]

    colors = ['#791a15','#f2b407','#181818']

    labels = ['2','3+','1']
    for j, df in enumerate([mut2, mut3plus, mut1]):
        plt.scatter(df['median_deltaG'],np.log2(df['median_norm_RNA:gDNA']),s=4,alpha=0.5,marker=markers[i],
                    label=msp.identifier,c=colors[j])
    
        out_of_bounds += len(df.loc[df['median_norm_RNA:gDNA']<y_limits_FL[0]])
        out_of_bounds += len(df.loc[df['median_norm_RNA:gDNA']>y_limits_FL[1]])
        
# plt.xlabel('deltaG')
# plt.ylabel('Log2 normalized RNA:DNA')
plt.ylim(y_limits_FL)
# plt.legend(loc='best')
plt.axvline(x=unmutated_dG,alpha=0.1,linestyle='--',color='k')
plt.axhline(y=0,alpha=0.1,linestyle='--',color='r')
plt.savefig('%s_thresh:%s_dG_vs_Ratio_collapsed_%s+_barcodes.pdf' % ('all_cdn_FL',str(msp_downstream_mutations.count_threshold), str(min_barcodes)),format='pdf')
plt.close()

records.write('Out of bounds in deltaG plot (full-length): %s\n' % str(out_of_bounds))