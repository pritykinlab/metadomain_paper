import numpy as np
import scipy
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
import matplotlib.pyplot as plt
import bbi
import pybedtools as pbt
import pickle


def make_chip_plot(full_clusters, peak_bedtool):
    n = len(full_clusters)
    delta = 2000
    fig, ax = plt.subplots()
    for u in np.unique(full_clusters):
        cluster_bedtool = pbt.BedTool([all_sig_regions[x] for x in range(n) if full_clusters[x]==u])
        chroms, ss, es = [], [],[]                
        for i in peak_bedtool.intersect(cluster_bedtool, u=True):
            chrom, s, e = i[:3]
            s, e = map(int, [s, e])
            chroms.append('chr' + chrom)
            s = (e+s)//2
            ss.append(s-delta)
            es.append(s+delta)
        if len(chroms) == 0:
            delta_plot = np.zeros(delta*2)
        else:
            tcon_vals = tcon_bw.stackup(chroms, ss, es)
            treg_vals = treg_bw.stackup(chroms, ss, es)
            delta_vals = treg_vals - tcon_vals
            delta_plot = np.nanmean(delta_vals, axis=0)
        ax.plot(delta_plot, label=f'{u}:{len(chroms)}')
    ax.set_title(f'{tf_name}:Treg - Tcon')
    ax.legend()
    return fig

peak_dict = {

    # 'GSE99319_p65' : ['GSE99319_p65_IP_in_Treg_COMBINED_peaks.narrowPeak', 'GSE99319_p65_IP_in_Tconv_COMBINED_peaks.narrowPeak'],
    # 'GSE69162_Satb1_' : ['GSE69162_Treg_Satb1_COMBINED_peaks.narrowPeak', 'GSE69162_Tcon_Satb1_COMBINED_peaks.narrowPeak'],
    # 'GSE40684_ets1' : ['GSE40684_treg_ets1_COMBINED_peaks.narrowPeak', 'GSE40684_cd4_ets1_COMBINED_peaks.narrowPeak', ],
    # 'GSE40684_cbfb' : ['GSE40684_treg_cbfb_COMBINED_peaks.narrowPeak', 'GSE40684_cd4_cbfb_COMBINED_peaks.narrowPeak'],
    # 'GSE40684_elf1' : ['GSE40684_treg_elf1_COMBINED_peaks.narrowPeak', 'GSE40684_cd4_elf1_COMBINED_peaks.narrowPeak'],
    # 'GSE40657_Foxo1_biot' : ['GSE40657_nTreg_Foxo1_Biotin_ChIPseq_COMBINED_peaks.narrowPeak', 'GSE40657_CD4+_Foxo1_Biotin_ChIPseq_COMBINED_peaks.narrowPeak'],
    # 'GSE40657_Foxo1_chip' : ['GSE40657_nTreg_Foxo1_Ab_ChIPseq_COMBINED_peaks.narrowPeak', 'GSE40657_CD4+_Foxo1_Ab_ChIPseq_COMBINED_peaks.narrowPeak'],
    'DRP003376_H3K27ac' : ['DRP003376_Treg_H3K27ac_COMBINED_peaks.narrowPeak', 'DRP003376_Tcon_H3K27ac_COMBINED_peaks.narrowPeak'],
    'DRP003376_SMC1a' : ['DRP003376_SMC1a_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_SMC1a_Tconv_COMBINED_peaks.narrowPeak'],
    'DRP003376_Satb1' : ['DRP003376_Satb1_TREG_COMBINED_peaks.narrowPeak', 'DRP003376_Satb1_TCON_COMBINED_peaks.narrowPeak'],
    'DRP003376_Runx1' : ['DRP003376_Runx1_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_Runx1_Tcon_COMBINED_peaks.narrowPeak'],
    'DRP003376_Mediator1' : ['DRP003376_Mediator1_COMBINED_peaks.narrowPeak', 'DRP003376_Mediator1_Tcon_COMBINED_peaks.narrowPeak'],
    'DRP003376_MBD' : ['DRP003376_MBD-seq_of_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_MBD-seq_of_Tconv_COMBINED_peaks.narrowPeak'],
    'DRP003376_H3K27me3' : ['DRP003376_H3K27me3_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_H3K27me3_Tconv_COMBINED_peaks.narrowPeak'],
    # 'DRP003376_H3K27ac' : ['DRP003376_H3K27ac_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_H3K27ac_Tcon_COMBINED_peaks.narrowPeak'],
    'DRP003376_H3K4me3' : ['DRP003376_H3K4me3_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_H3K4me3_Tcon_COMBINED_peaks.narrowPeak'],
    'DRP003376_H3K4me1' : ['DRP003376_H3K4me1_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_H3K4me1_Tcon_COMBINED_peaks.narrowPeak'],
    'DRP003376_Ets1' : ['DRP003376_Ets1_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_Ets1_Tcon_COMBINED_peaks.narrowPeak'],
    'DRP003376_CREB' : ['DRP003376_CREB_Treg_COMBINED_peaks.narrowPeak', 'DRP003376_CREB_Tconv_COMBINED_peaks.narrowPeak'],
    'DRP003376_Bcl11b' : ['DRP003376_Bcl11b_COMBINED_peaks.narrowPeak', 'DRP003376_Bcl11b_Tcon_COMBINED_peaks.narrowPeak'],
}

bw_pref = '../bws/'
all_sig_regions = pickle.load( open( "all_sig_regions.p", "rb" ) )
full_clusters = pickle.load( open( "full_clusters.p", "rb" ) )
pref = '../narrowPeaks/'

for tf_name, vals in peak_dict.items():
    treg_peak_file, tcon_peak_file = vals

    treg_bedtool = pbt.BedTool(pref + treg_peak_file)
    tcon_bedtool = pbt.BedTool(pref + tcon_peak_file)
    
    both = []
    tcon_without_treg = tcon_bedtool.subtract(treg_bedtool, A=True)
    for i in tcon_without_treg:
        both.append(tuple(i[:3]))
    treg_without_tcon = treg_bedtool.subtract(tcon_bedtool, A=True)
    for i in treg_without_tcon:
        both.append(tuple(i[:3]))
    extreme_bedtool = pbt.BedTool(both).sort().merge().saveas()

    both = []
    tcon_without_treg = tcon_bedtool
    for i in tcon_without_treg:
        both.append(tuple(i[:3]))
    treg_without_tcon = treg_bedtool
    for i in treg_without_tcon:
        both.append(tuple(i[:3]))
    union_bedtool = pbt.BedTool(both).sort().merge().saveas()

    both = []
    tcon_without_treg = tcon_bedtool.intersect(treg_bedtool, u=True)
    for i in tcon_without_treg:
        both.append(tuple(i[:3]))
    treg_without_tcon = treg_bedtool.intersect(tcon_bedtool, u=True)
    for i in treg_without_tcon:
        both.append(tuple(i[:3]))
    intersection_bedtool = pbt.BedTool(both).sort().merge().saveas()
    
    print(tcon_peak_file, treg_peak_file)
    tcon_base = tcon_peak_file.split('_peaks')[0]
    treg_base = treg_peak_file.split('_peaks')[0]
    print(tcon_base, treg_base)
    tcon_bw = bbi.open(bw_pref + tcon_base + ".bw")
    treg_bw = bbi.open(bw_pref + treg_base + ".bw")

    pbt.BedTool(all_sig_regions).intersect(extreme_bedtool, c=True)
    pbt.BedTool(all_sig_regions).intersect(intersection_bedtool, c=True)
    pbt.BedTool(all_sig_regions).intersect(union_bedtool, c=True)

    fig = make_chip_plot(full_clusters, extreme_bedtool)
    fig.savefig(f'./plots/{tf_name}_extreme.jpg')
    plt.close(fig)

    fig = make_chip_plot(full_clusters, intersection_bedtool)
    fig.savefig(f'./plots/{tf_name}_intersection.jpg')
    plt.close(fig)

    fig = make_chip_plot(full_clusters, union_bedtool)
    fig.savefig(f'./plots/{tf_name}_union.jpg')
    plt.close(fig)

