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


bw_dict = {
    'Ep300ChIP_Tconv' : ['DRA010814_Ep300ChIP_Treg_periphery_WT_COMBINED.bw', 'DRA010814_Ep300ChIP_Tconv_periphery_WT_COMBINED.bw',],
    'Stat5ChIP_Tconv' : ['DRA010814_Stat5ChIP_Treg_periphery_WT_COMBINED.bw', 'DRA010814_Stat5ChIP_Tconv_periphery_WT_COMBINED.bw',],
    'Tet2ChIP_Tconv' : ['DRA010814_Tet2ChIP_Treg_periphery_WT_COMBINED.bw', 'DRA010814_Tet2ChIP_Tconv_periphery_WT_COMBINED.bw',],
    'cJun' : ['GSE154680_cJun-CUTRUN_resting_Treg_COMBINED.bw', 'GSE154680_cJun-CUTRUN_resting_Tcon_COMBINED.bw',],
    'IRF4' : ['GSE154680_IRF4-CUTRUN_resting_Treg_COMBINED.bw', 'GSE154680_IRF4-CUTRUN_resting_Tcon_COMBINED.bw',],
    'Lef1' : ['GSE154680_Lef1-CUTRUN_resting_Treg_COMBINED.bw', 'GSE154680_Lef1-CUTRUN_resting_Tcon_COMBINED.bw',],
    'TCF1' : ['GSE154680_TCF1-CUTRUN_resting_Treg_COMBINED.bw', 'GSE154680_TCF1-CUTRUN_resting_Tcon_COMBINED.bw',],

    'Bcl11b_Tcon' : ['DRP003376_Bcl11b_COMBINED.bw', 'DRP003376_Bcl11b_Tcon_COMBINED.bw'],
    'CREB_Tconv' : ['DRP003376_CREB_Treg_COMBINED.bw', 'DRP003376_CREB_Tconv_COMBINED.bw'],
    'Ets1_Tcon' : ['DRP003376_Ets1_Treg_COMBINED.bw', 'DRP003376_Ets1_Tcon_COMBINED.bw'],
    'MBD-seq_of_Tconv' : ['DRP003376_MBD-seq_of_Treg_COMBINED.bw',  'DRP003376_MBD-seq_of_Tconv_COMBINED.bw'],
    'Mediator1_Tcon' : ['DRP003376_Mediator1_COMBINED.bw', 'DRP003376_Mediator1_Tcon_COMBINED.bw'],
    'Runx1_Tcon' : ['DRP003376_Runx1_Treg_COMBINED.bw', 'DRP003376_Runx1_Tcon_COMBINED.bw'],
    'Satb1_TCON' : ['DRP003376_Satb1_TREG_COMBINED.bw', 'DRP003376_Satb1_TCON_COMBINED.bw'],
    'SMC1a_Tconv' : ['DRP003376_SMC1a_Treg_COMBINED.bw', 'DRP003376_SMC1a_Tconv_COMBINED.bw'],
    'CD4+_Foxo1_Ab_ChIPseq' : ['GSE40657_nTreg_Foxo1_Ab_ChIPseq_COMBINED.bw', 'GSE40657_CD4+_Foxo1_Ab_ChIPseq_COMBINED.bw'],
    'CD4+_Foxo1_Biotin_ChIPseq' : ['GSE40657_nTreg_Foxo1_Biotin_ChIPseq_COMBINED.bw', 'GSE40657_CD4+_Foxo1_Biotin_ChIPseq_COMBINED.bw'],
    'cd4_elf1' : ['GSE40684_treg_elf1_COMBINED.bw', 'GSE40684_cd4_elf1_COMBINED.bw'],
    'cd4_cbfb' : ['GSE40684_treg_cbfb_COMBINED.bw', 'GSE40684_cd4_cbfb_COMBINED.bw'],
    'cd4_ets1' : ['GSE40684_treg_ets1_COMBINED.bw', 'GSE40684_cd4_ets1_COMBINED.bw'],
    'p65_IP_in_Tconv' : ['GSE99319_p65_IP_in_Treg_COMBINED.bw', 'GSE99319_p65_IP_in_Tconv_COMBINED.bw'],
    'tn_h3k27me3' : ['GSE55753_treg_h3k27me3_COMBINED.bw', 'GSE55753_tn_h3k27me3_COMBINED.bw',],
    'H3K27me3_Tconv' : ['DRP003376_H3K27me3_Treg_COMBINED.bw', 'DRP003376_H3K27me3_Tconv_COMBINED.bw',],
    'H3K27me3_immature_CD4SP' : ['DRP003376_H3K27me3_Treg_COMBINED.bw', 'DRP003376_H3K27me3_ChIP-seq_of_immature_CD4SP_COMBINED.bw',],
    'H3K27me3_tTreg_precursor' : ['DRP003376_H3K27me3_Treg_COMBINED.bw', 'DRP003376_H3K27me3_ChIP-seq_of_tTreg_precursor_COMBINED.bw',],
    'H3K27me3_tTreg' : ['DRP003376_H3K27me3_Treg_COMBINED.bw', 'DRP003376_H3K27me3_ChIP-seq_of_tTreg_COMBINED.bw',],
    'H3K4me1_Tcon' : ['DRP003376_H3K4me1_Treg_COMBINED.bw', 'DRP003376_H3K4me1_Tcon_COMBINED.bw',],
    'H3K4me1_Cd4CreSatb1CKO_immature_CD4SP' : ['DRP003376_H3K4me1_Treg_COMBINED.bw', 'DRP003376_H3K4me1_ChIP-seq_of_Cd4CreSatb1CKO_immature_CD4SP_COMBINED.bw',],
    'H3K4me1_tTreg_precursor' : ['DRP003376_H3K4me1_Treg_COMBINED.bw', 'DRP003376_H3K4me1_ChIP-seq_of_tTreg_precursor_COMBINED.bw',],
    'H3K4me1_tTreg' : ['DRP003376_H3K4me1_Treg_COMBINED.bw', 'DRP003376_H3K4me1_ChIP-seq_of_tTreg_COMBINED.bw',],
    'H3K4me3_Tcon' : ['DRP003376_H3K4me3_Treg_COMBINED.bw', 'DRP003376_H3K4me3_Tcon_COMBINED.bw',],
    'H3K4me3_immature_CD4SP' : ['DRP003376_H3K4me3_Treg_COMBINED.bw', 'DRP003376_H3K4me3_ChIP-seq_of_immature_CD4SP_COMBINED.bw',],
    'H3K4me3_tTreg_precursor' : ['DRP003376_H3K4me3_Treg_COMBINED.bw', 'DRP003376_H3K4me3_ChIP-seq_of_tTreg_precursor_COMBINED.bw',],
    'H3K4me3_tTreg' : ['DRP003376_H3K4me3_Treg_COMBINED.bw', 'DRP003376_H3K4me3_ChIP-seq_of_tTreg_COMBINED.bw',],
    'Satb1_TCON' : ['DRP003376_Satb1_TREG_COMBINED.bw', 'DRP003376_Satb1_TCON_COMBINED.bw',],
    'Satb1_CD4SP' : ['DRP003376_Satb1_TREG_COMBINED.bw', 'DRP003376_Satb1_ChIP-seq_of_CD4SP_COMBINED.bw',],
    'H3K27ac_Tcon' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_Tcon_COMBINED.bw',],
    'Satb1_tTreg_precurs_tTreg' : ['DRP003376_Satb1_TREG_COMBINED.bw', 'DRP003376_Satb1_ChIP-seq_of_tTreg_precursor_and_tTreg_COMBINED.bw',],
    'H3K27ac_activated_Tconv' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_activated_Tconv_COMBINED.bw',],
    'H3K27ac_activated_Treg' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_activated_Treg_COMBINED.bw',],
    'H3K27ac_Cd4CreSatb1CKO_immature_CD4SP' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_Cd4CreSatb1CKO_immature_CD4SP_COMBINED.bw',],
    'H3K27ac_Cd4CreSatb1CKO_tTreg_precursor' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_Cd4CreSatb1CKO_tTreg_precursor_COMBINED.bw',],
    'H3K27ac_immature_CD4SP' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_immature_CD4SP_COMBINED.bw',],
    'H3K27ac_tTreg_precursor' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_tTreg_precursor_COMBINED.bw',],
    'H3K27ac_tTreg' : ['DRP003376_H3K27ac_ChIP-seq_of_Treg_COMBINED.bw', 'DRP003376_H3K27ac_ChIP-seq_of_tTreg_COMBINED.bw',]
}

def make_order2(matrix):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method='average', metric='cosine')
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)

        order = dendro['leaves']
    else:
        order = np.arange(len(matrix))
    # order = dendro['leaves']
    return order        

def make_delta_order(peak_dict, gene_bedtool):
    n = len(peak_dict)
    delta = 1000
    order_dict = {}
    for u in peak_dict:
        chroms, ss, es = [], [],[]                
        print("Doing", u)
        for i in peak_dict[u]:
            chrom, s, e = i[:3]
            s, e = map(int, [s, e])
            chroms.append('chr' + chrom)
            ss.append(s-5000)
            es.append(e+5000)
        if len(chroms) <= 10:
            delta_plot = np.zeros(delta*2)
        else:
            cond_vals = cond_bw.stackup(chroms, ss, es, bins=delta*2)
            treg_vals = treg_bw.stackup(chroms, ss, es, bins=delta*2)
            delta_vals = treg_vals - cond_vals
            delta_vals[np.isnan(delta_vals)] = 0
            delta_vals[np.isinf(delta_vals)] = 0
            delta_vals[:, 0] += 1e-10
            delta_order = make_order2(delta_vals)
            order_dict[u] = delta_order
    return order_dict


def make_loop_chip_plot(peak_dict, gene_bedtool):
    n = len(peak_dict)
    delta = 1000
    delta_dict = {}
    cond_dict = {}    
    treg_dict = {}
    for u in peak_dict:
        chroms, ss, es = [], [],[]                
        for i in peak_dict[u]:
            chrom, s, e = i[:3]
            s, e = map(int, [s, e])
            chroms.append('chr' + chrom)
            ss.append(s-5000)
            es.append(e+5000)
        if len(chroms) <= 10:
            delta_plot = np.zeros(delta*2)
        else:
            cond_vals = cond_bw.stackup(chroms, ss, es, bins=delta*2)
            treg_vals = treg_bw.stackup(chroms, ss, es, bins=delta*2)

            delta_vals = treg_vals - cond_vals
            delta_dict[u] = delta_vals
            treg_dict[u] = treg_vals
            cond_dict[u] = cond_vals
    return treg_dict, cond_dict, delta_dict


bw_pref = '../bws/'
all_tcon_v_treg = pbt.BedTool('Tcon_all.csv')
peak_file = all_tcon_v_treg

tcon_loops = pbt.BedTool("tcon_ancs.csv")
treg_loops = pbt.BedTool("treg_ancs.csv")
nonsig_loops = pbt.BedTool("nonsignificant_ancs.csv")



loop_dict = {
    'tcon_loops' : all_tcon_v_treg.intersect(tcon_loops, u=True),
    'treg_loops' : all_tcon_v_treg.intersect(treg_loops, u=True),
    'nonsig_loops' : all_tcon_v_treg.intersect(nonsig_loops, u=True),
}


n_clusts = len(loop_dict)
n_bws = len(bw_dict)
slopsize = 50_000

for c, (tf_name, vals) in enumerate(bw_dict.items()):
    treg_peak_file, cond_peak_file = vals
    cond_bw = bbi.open(bw_pref + cond_peak_file)
    treg_bw = bbi.open(bw_pref + treg_peak_file)

del cond_bw
del treg_bw


import time

treg_bw = bbi.open(bw_pref + 'DRP003376_Treg_H3K27ac_COMBINED.bw')
cond_bw = bbi.open(bw_pref + 'DRP003376_Tcon_H3K27ac_COMBINED.bw')
delta_order_dict = make_delta_order(loop_dict, peak_file)



line_fig, line_ax = plt.subplots(n_bws, 3, figsize=(3*4, 4*n_bws))

for c, (tf_name, vals) in enumerate(bw_dict.items()):
    delta_fig, delta_ax = plt.subplots(1, n_clusts, figsize=(4*n_clusts, 4))
    treg_fig, raw_treg_ax = plt.subplots(1, n_clusts, figsize=(4*n_clusts, 4))
    condition_fig, raw_condition_ax = plt.subplots(1, n_clusts, figsize=(4*n_clusts, 4))

    treg_peak_file, cond_peak_file = vals

    cond_bw = bbi.open(bw_pref + cond_peak_file)
    treg_bw = bbi.open(bw_pref + treg_peak_file)

    t1 = time.time()
    treg_dict, cond_dict, delta_dict = make_loop_chip_plot(loop_dict, peak_file)
    t1f = time.time()
    print("Make dicts", t1f-t1)
    for u, i in enumerate(loop_dict):
        delta_order = delta_order_dict[i]        

        delta_ax[u].matshow(delta_dict[i][delta_order, :], aspect = 'auto', cmap=cm.bwr, vmin=-2, vmax=2)
        raw_treg_ax[u].matshow(treg_dict[i][delta_order, :], aspect = 'auto', cmap=cm.bwr, vmin=0, vmax=10)
        raw_condition_ax[u].matshow(cond_dict[i][delta_order, :], aspect = 'auto', cmap=cm.bwr, vmin=0, vmax=10)

        delta_ax[u].set_title(f'delta {tf_name}, cluster={i}')
        raw_treg_ax[u].set_title(f'treg {tf_name}, cluster={i}')
        raw_condition_ax[u].set_title(f'cond {tf_name}, cluster={i}')

        line_ax[c, 0].plot(np.nanmean(delta_dict[i], axis=0), label=f'{i}')
        line_ax[c, 1].plot(np.nanmean(treg_dict[i], axis=0), label=f'{i}')
        line_ax[c, 2].plot(np.nanmean(cond_dict[i], axis=0), label=f'{i}')
        
    line_ax[c, 0].legend()
    line_ax[c, 0].set_title(f"Delta {tf_name}")
    line_ax[c, 1].set_title(f"Treg {tf_name}")
    line_ax[c, 2].set_title(f"Cond {tf_name}")
    delta_fig.tight_layout()
    treg_fig.tight_layout()
    condition_fig.tight_layout()
    delta_fig.savefig(f'./plots/delta_{tf_name}.jpg')
    treg_fig.savefig(f'./plots/treg_{tf_name}.jpg')
    condition_fig.savefig(f'./plots/tcon_{tf_name}.jpg')

    plt.close(delta_fig)
    plt.close(treg_fig)
    plt.close(condition_fig)


line_fig.tight_layout()

line_fig.savefig('./plots/lines.jpg')




