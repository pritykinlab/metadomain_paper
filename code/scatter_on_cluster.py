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
from matplotlib.markers import MarkerStyle

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

def differential_hic_scatterplot(all_loops, cool_treg, cool_tconv):
	tconv_scores = []
	treg_scores = []

	wsz = 4
	for i in all_loops:
		l1, l2 = i[:3], i[3:6]
		
		l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
		l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)
		
		val_tconv = cool_tconv.matrix(balance=True).fetch(l1, l2).sum()
		val_treg = cool_treg.matrix(balance=True).fetch(l1, l2).sum()
		tconv_scores.append(val_tconv)
		treg_scores.append(val_treg)
	tconv_scores, treg_scores = np.asarray(tconv_scores), np.asarray(treg_scores)
	return tconv_scores, treg_scores

all_ancs = pbt.BedTool("full_ancs.csv")
all_loops = pbt.BedTool("full_loops.loops")
all_tcon_v_treg = pbt.BedTool('Tcon_all.csv')

import cooler
cool_treg = cooler.Cooler(f'./Treg_all.mcool::resolutions/5000')
cool_tconv = cooler.Cooler(f'./Tconv_all.mcool::resolutions/5000')

# tcon_scores, treg_scores = differential_hic_scatterplot(all_loops, cool_treg, cool_tconv)
import pickle
# pickle.dump(tcon_scores, open( "all_tcon_scores.p", "wb" ) )
# pickle.dump(treg_scores, open( "all_treg_scores.p", "wb" ) )

tcon_scores = pickle.load( open( "all_tcon_scores.p", "rb" ) )
treg_scores = pickle.load( open( "all_treg_scores.p", "rb" ) )

ancs = []
chroms, ss, es = [], [],[]         
delta = 500       
for i in all_ancs.intersect(all_tcon_v_treg, wo=True):
	anc = tuple(i[:3])
	atac_peak = tuple(i[3:6])
	chrom, s, e = atac_peak[:3]
	s, e = map(int, [s, e])
	chroms.append('chr' + chrom)
	ss.append(s)
	es.append(e)
	ancs.append(anc)

bw_pref = '../bws/'
n = len(bw_dict)
fig, ax = plt.subplots(n//4, 4, figsize=(4*4, 4*n//4))
ax = np.ravel(ax)
for c, (tf_name, vals) in enumerate(bw_dict.items()):
	treg_peak_file, cond_peak_file = vals

	cond_bw = bbi.open(bw_pref + cond_peak_file)
	treg_bw = bbi.open(bw_pref + treg_peak_file)

	cond_vals = np.nanmean(cond_bw.stackup(chroms, ss, es, bins=delta*2), axis=1)
	treg_vals = np.nanmean(treg_bw.stackup(chroms, ss, es, bins=delta*2), axis=1)

	anc_to_vals = {} 
	for c2, i in enumerate(ancs):
		anc_to_vals.setdefault(i, [])
		anc_to_vals[i].append(cond_vals[c2] - treg_vals[c2])
	for i in all_ancs:
		anc = tuple(i[:3])
		if not anc_to_vals.get(anc):
			anc_to_vals[anc] = 0 
	l_vals = []
	r_vals = []		
	for i in all_loops:
		l1, l2 = tuple(i[:3]), tuple(i[3:6])
		L_val = np.nanmean(anc_to_vals[l1])
		R_val = np.nanmean(anc_to_vals[l2])
		l_vals.append(L_val)
		r_vals.append(R_val)

	a = ax[c]

	a.scatter(tcon_scores, treg_scores, c = l_vals,
				marker=MarkerStyle('o', fillstyle='left'),
				cmap=cm.bwr, vmin=-2, vmax=2, label='left')

	a.scatter(tcon_scores, treg_scores, c = r_vals,
				marker=MarkerStyle('o', fillstyle='right'),
				cmap=cm.bwr, vmin=-2, vmax=2, label='right')

	a.set_title(tf_name)
fig.savefig('./plots/scatter.jpg')
