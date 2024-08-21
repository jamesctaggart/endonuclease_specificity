import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pdb
from mutscan.sequtils import calculate_GC
from collections import defaultdict

def get_fig_axes(axes=None):
    """Retrieve figure and axes from `axes`. If `axes` is None, both.

    Used as a helper function for replotting atop existing axes, by functions
    defined in :mod:`plastid.plotting.plots`.

	From plastid 0.4.8 (https://plastid.readthedocs.io/en/latest/)

    Parameters
    ----------
    axes : :class:`matplotlib.axes.Axes` or `None`
        Axes in which to place plot. If `None`, a new figure is generated.


    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Parent figure of axes

    :class:`matplotlib.axes.Axes`
        Axes containing plot
    """
    if axes is None:
        fig = plt.figure()
        ax = plt.gca()
    else:
        ax = axes
        fig = ax.figure

    return fig, ax

def cdf(data,axes=None,**kwargs):
    fig, ax = get_fig_axes(axes)
    # NOTE: This function may produce weird results if bin size is too small
    num_bins = 1000000
    counts, bin_edges = np.histogram(data, bins=num_bins, density=True)
    cdf = np.cumsum(counts)
    ax.plot(bin_edges[1:], cdf/cdf[-1], **kwargs)


def snp_boxplot(SNP_ratios, ylim, figname, figsize, topbar=True, vert_stripe = False, hlines = [],
            logy=False, show_wt=True, N_y=0, zero_count=True, record_zeros=None, **kwargs):
    
    boxprops = dict(linestyle='-', linewidth=1, color='black')
    medianprops = dict(linestyle='-', linewidth=3, color='k')
    whiskerprops = dict(linestyle='--', color='black')


    if not topbar:
        fig, ax = plt.subplots(figsize=figsize)

        if not logy:
            plt.axhspan(np.percentile(SNP_ratios['wt'],75),
                        np.percentile(SNP_ratios['wt'],25), facecolor='gray', alpha=0.25)

        else:
            plt.axhspan(np.log2(np.percentile(SNP_ratios['wt'],75)),
                np.log2(np.percentile(SNP_ratios['wt'],25)), facecolor='gray', alpha=0.25)

        if not show_wt:
            SNP_ratios.pop('wt')

        if not logy:
            bplot = ax.boxplot(SNP_ratios.values(), labels=SNP_ratios.keys(), whis=[5, 95],
            showfliers=False,medianprops=medianprops,boxprops=boxprops,whiskerprops=whiskerprops,patch_artist=True,zorder=3, **kwargs)
            for ii, x in enumerate(SNP_ratios.values()):
                ax.text(ii + 0.8, 2.1, str(len(x)), rotation=90)
        else:
            bplot = ax.boxplot([np.log2(x) for x in SNP_ratios.values()], labels=SNP_ratios.keys(), whis=[5, 95],
            showfliers=False,medianprops=medianprops,boxprops=boxprops,whiskerprops=whiskerprops, patch_artist=True,zorder=3, **kwargs)
            for ii, x in enumerate(SNP_ratios.values()):
                ax.text(ii + 0.8, N_y, str(len(x)), rotation=90)

        ax.tick_params(axis='x', labelrotation=90)

    else:
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, hspace=0, height_ratios=[1, 6])
        axs = gs.subplots(sharex=True, sharey=False)

        if not logy:
            plt.axhspan(np.percentile(SNP_ratios['wt'],75),
                        np.percentile(SNP_ratios['wt'],25), facecolor='gray', alpha=0.25)

        else:
            plt.axhspan(np.log2(np.percentile(SNP_ratios['wt'],75)),
                np.log2(np.percentile(SNP_ratios['wt'],25)), facecolor='gray', alpha=0.25)

        if not show_wt:
            SNP_ratios.pop('wt')

        if not logy:
            bplot = axs[1].boxplot(SNP_ratios.values(), labels=SNP_ratios.keys(), whis=[5, 95],
            showfliers=False,medianprops=medianprops,boxprops=boxprops,whiskerprops=whiskerprops,patch_artist=True,zorder=3, **kwargs)

        else:
            bplot = axs[1].boxplot([np.log2(x) for x in SNP_ratios.values()], labels=SNP_ratios.keys(), whis=[5, 95],
            showfliers=False,medianprops=medianprops,boxprops=boxprops,whiskerprops=whiskerprops, patch_artist=True,zorder=3, **kwargs)

        bar_colors = []  # Define bar colors that are grouped by threes 
        for ii in range(0, len(SNP_ratios.values())):
            if (np.floor(ii / 3.)) % 2 == 0:
                bar_colors.append('DarkGray')
            else:
                bar_colors.append('LightGray')

        axs[0].bar(range(1, len(SNP_ratios.values())+1),[len(x) for x in SNP_ratios.values()],color=bar_colors)

        axs[1].tick_params(axis='x', labelrotation=90)
        axs[0].spines[['right', 'top']].set_visible(False)

    if vert_stripe:
        for i in np.arange(0,(len(SNP_ratios.values())-show_wt)/3,2): # the -show_wt means that if WT is shown you don't draw an extra box for it
            plt.axvspan(0.5+3*i, 0.5+3*(i+1), facecolor='gray', alpha=0.05,zorder=2)
    if len(hlines)>0:
        for y in hlines:
            plt.axhline(y=y, color='k', linestyle='--', alpha=0.2,zorder=2)

    for ii, patch in enumerate(bplot['boxes']):
        if (np.floor(ii / 3.)) % 2 == 0:
            patch.set_facecolor('LightSeaGreen')
        else:
            patch.set_facecolor('PaleTurquoise')

    if record_zeros is not None:
        y_offset = 0.1
        zero_val_count = {k:len([x for x in v if x==0.0]) for k,v in SNP_ratios.items()}  # Count the number of zero values in the SNP_ratios
        for k,v in {k:v for k,v in zero_val_count.items() if v!=0}.items():
            record_zeros.write(k+'\t'+str(v)+'\n')
            
    # plt.ylim(0, 2)
    plt.ylabel('RNA / DNA Ratio')
    plt.ylim(ylim)
    plt.savefig(figname,format='pdf')
    plt.close()


def snp_boxplot_topbars(SNP_ratios,figname, figsize, vert_stripe = False, hlines = [],
            logy=False,show_wt=True, **kwargs):

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, hspace=0, height_ratios=[1, 6])
    axs = gs.subplots(sharex=True, sharey=False)

    boxprops = dict(linestyle='-', linewidth=1, color='black')
    medianprops = dict(linestyle='-', linewidth=3, color='k')
    whiskerprops = dict(linestyle='--', color='black')

    # breakpoint()
    if not logy:
        plt.axhspan(np.percentile(SNP_ratios['wt'],75),
                    np.percentile(SNP_ratios['wt'],25), facecolor='gray', alpha=0.25)

    else:
        plt.axhspan(np.log2(np.percentile(SNP_ratios['wt'],75)),
            np.log2(np.percentile(SNP_ratios['wt'],25)), facecolor='gray', alpha=0.25)

    if not show_wt:
        SNP_ratios.pop('wt')

    if not logy:
        bplot = axs[1].boxplot(SNP_ratios.values(), labels=SNP_ratios.keys(), whis=[5, 95],
        showfliers=False,medianprops=medianprops,boxprops=boxprops,whiskerprops=whiskerprops,patch_artist=True,zorder=3, **kwargs)

    else:
        bplot = axs[1].boxplot([np.log2(x) for x in SNP_ratios.values()], labels=SNP_ratios.keys(), whis=[5, 95],
        showfliers=False,medianprops=medianprops,boxprops=boxprops,whiskerprops=whiskerprops, patch_artist=True,zorder=3, **kwargs)

    axs[0].bar(range(1, len(SNP_ratios.values())+1),[len(x) for x in SNP_ratios.values()],color='gray')

    if vert_stripe:
        for i in np.arange(0,(len(SNP_ratios.values())-show_wt)/3,2): # the -show_wt means that if WT is shown you don't draw an extra box for it
            plt.axvspan(0.5+3*i, 0.5+3*(i+1), facecolor='gray', alpha=0.05,zorder=2)
    if len(hlines)>0:
        for y in hlines:
            plt.axhline(y=y, color='k', linestyle='--', alpha=0.2,zorder=2)

    for ii, patch in enumerate(bplot['boxes']):
        if (np.floor(ii / 3.)) % 2 == 0:
            patch.set_facecolor('LightSeaGreen')
        else:
            patch.set_facecolor('PaleTurquoise')

    axs[1].tick_params(axis='x', labelrotation=90)

    # plt.ylim(0, 2)
    # plt.ylabel('RNA / DNA Ratio')
    plt.savefig(figname,format='pdf')
    plt.close()


def double_mutant_heatmap(doublemut_df,axes=None, zero_position=0,**kwargs):
    fig, ax = get_fig_axes(axes)

    mask = doublemut_df.isnull()
    sns.heatmap(doublemut_df, mask=mask, vmin=-5, vmax=2,
                cmap='RdBu')  # For color see: https://www.r-graph-gallery.com/38-rcolorbrewers-palettes
    plt.savefig('doublemut_cleanbg_headmap.pdf', format='pdf', dpi=300)
    plt.close()

def make_barcode_bias_plots(barcode_bias_file, dataset_id, out_base):

    # Load in barcode bias file
    biases = {}
    with open(barcode_bias_file,'r') as f:
        for line in f:
            fields = line[:-1].split()
            biases[fields[0]] = float(fields[1])

    # Make plot of GC content vs barcode bias
    fig, ax = plt.subplots()
    make_bias_vs_GC(biases, axes=ax)
    plt.savefig('%s/%s_bias_vs_GC' % (out_base, dataset_id))
    plt.close()

    make_bias_vs_nucleotide_content(biases, dataset_id, out_base)


def make_bias_vs_GC(bc_dict, axes=None):
    fig, ax = get_fig_axes(None)

    binned = defaultdict(list)
    bc_dict_values = list(bc_dict.values())
    bc_dict_keys = list(bc_dict.keys())

    for pair in sorted(zip(bc_dict_values, bc_dict_keys)):
        binned[calculate_GC(pair[1])].append(pair[0])
            
    sorted_binned = sorted(binned.items())
    labels = [str(x)[:4] for x, _ in sorted_binned]
    data = [y for _, y in sorted_binned]

    ax.boxplot(data,whis=[5,95],labels=labels)

    ax.set_xlabel('Fraction Gor C')
    plt.xticks(rotation='vertical')
    ax.set_ylim([0,2])
    ax.set_label('Normalized RNA/DNA')
    # plt.savefig('wt_barcode_ratio_vs_GC.pdf',format='pdf')

def make_bias_vs_nucleotide_content(bc_dict,dataset_id,out_base,axes=None):
    
    fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

    fig.tight_layout(pad=3.0)
    bc_dict_values = list(bc_dict.values())
    bc_dict_keys = list(bc_dict.keys())

    binnedA = defaultdict(list)
    for pair in sorted(zip(bc_dict_values, bc_dict_keys)):
        binnedA[pair[1].count('A')].append(pair[0])
    binnedT = defaultdict(list)
    for pair in sorted(zip(bc_dict_values, bc_dict_keys)):
        binnedT[pair[1].count('T')].append(pair[0])
    binnedC = defaultdict(list)
    for pair in sorted(zip(bc_dict_values, bc_dict_keys)):
        binnedC[pair[1].count('C')].append(pair[0])
    binnedG = defaultdict(list)
    for pair in sorted(zip(bc_dict_values, bc_dict_keys)):
        binnedG[pair[1].count('G')].append(pair[0])
        
    plt.xticks(rotation='vertical')
    datasets = (binnedA,binnedT,binnedC,binnedG)
    nt = ['A','T','C','G']

    for ii,ax in enumerate((ax1,ax2,ax3,ax4)):
        ds = datasets[ii]
                
        sorted_binned = sorted(ds.items())
        labels = [str(x)[:4] for x, _ in sorted_binned]
        data = [y for _, y in sorted_binned]

        ax.boxplot(data,whis=[5,95],labels=labels)

        ax.set_xlabel('count %s' % nt[ii])
        ax.set_ylim([0,2])
        ax.set_ylabel('Normalized RNA/DNA')

    plt.savefig('%s/%s_bias_vs_sequence_content' % (out_base, dataset_id))
    plt.close()