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
import scipy
import scipy.ndimage
arr = np.asarray

def subslice(sl, d):
    s = slice(sl.start-d, sl.stop-d)
    return s


def aggregate_slices(cool_dict, exp_dict, all_ind_to_region, all_connections,  chrom_to_start, chrom_to_end, chromsizes, d, dist_co=100):
    exp_vals = { x: [] for x in cool_dict }
    bias_vals = { x: [] for x in cool_dict }
    m_vals = { x: [] for x in cool_dict }
    inds = []

    labels = scipy.ndimage.label(all_connections)[0]
    slices = scipy.ndimage.find_objects(labels)

    differential_treg_tcon = []
    for chrom in chromsizes:    
        if ('M' in chrom) or ('Y' in chrom) or ('X' in chrom):
            continue
        for c, (cond, cool) in enumerate(cool_dict.items()):
            m = cool.matrix(balance=False).fetch(chrom)
            biasvec = cool.bins().fetch(chrom)['weight']
            biasmat = np.outer(biasvec, biasvec)

            rows, cols = np.indices(m.shape)
            indmat = cols-rows
            goodinds = indmat > 2
            
            exp_df = exp_dict[cond]
            exp_mat = get_exp_from_indmat(exp_df, chrom, indmat, goodinds)
            for sl in slices:
                sl1, sl2 = sl
                
                if abs(sl1.start - sl2.start > dist_co) or (sl1.start - sl2.start < 0):
                    continue
                chrom2, _, _ = all_ind_to_region[sl1.start]
                d = chrom_to_start[chrom2]
                sl1_new, sl2_new = subslice(sl1, d), subslice(sl2, d)
                
                if chrom == chrom2:
                    exp = np.sum((exp_mat[sl1_new, sl2_new])*(all_connections[sl1, sl2]))
                    bias = np.sum((biasmat[sl1_new, sl2_new])*(all_connections[sl1, sl2]))
                    val = np.nansum((m[sl1_new, sl2_new])*(all_connections[sl1, sl2]))
                    
                    exp_vals[cond].append(exp)
                    bias_vals[cond].append(bias)
                    m_vals[cond].append(val)
                    if c == 0:
                        inds.append(f'{sl1.start}_{sl1.stop}_{sl2.start}_{sl2.stop}')
                else:
                    pass
    inds = arr(inds)
    for cond in exp_vals:
        exp_vals[cond] = arr(exp_vals[cond])
    for cond in bias_vals:
        bias_vals[cond] = arr(bias_vals[cond])
    for cond in m_vals:
        m_vals[cond] = arr(m_vals[cond])
    return m_vals, exp_vals, bias_vals, inds, differential_treg_tcon

def aggregate_reads(cool_dict, exp_dict, all_ind_to_region,  chrom_to_start, chrom_to_end, chromsizes, places, deseq_pval_mat_2, d):
    exp_vals = { x: [] for x in cool_dict }
    bias_vals = { x: [] for x in cool_dict }
    exp_vals2 = { x: [] for x in cool_dict }
    bias_vals2 = { x: [] for x in cool_dict }

    m_vals = { x: [] for x in cool_dict }
    inds = []
    differential_treg_tcon = []
    for chrom in chromsizes:    
        if ('M' in chrom) or ('Y' in chrom) or ('X' in chrom):
            continue
        for c, (cond, cool) in enumerate(cool_dict.items()):
            m = cool.matrix(balance=False).fetch(chrom)
            biasvec = cool.bins().fetch(chrom)['weight']
            biasmat = np.outer(biasvec, biasvec)

            rows, cols = np.indices(m.shape)
            indmat = cols-rows
            goodinds = indmat > 2
            
            exp_df = exp_dict[cond]
            exp_mat = get_exp_from_indmat(exp_df, chrom, indmat, goodinds)
            for place in places:
                i1, i2 = place
                if (i2-i1 > 100) or (i2-i1) < 0:
                    continue
                chrom2, _, _ = all_ind_to_region[i1]
                ind1, ind2 = i1-chrom_to_start[chrom2], i2-chrom_to_start[chrom2]
                if chrom == chrom2:
                    exp = np.sum(exp_mat[ind1-d:ind1+d+1, ind2-d:ind2+d+1])
                    bias = np.sum(biasmat[ind1-d:ind1+d+1, ind2-d:ind2+d+1])

                    val = np.nansum(m[ind1-d:ind1+d+1, ind2-d:ind2+d+1])
                    
                    exp_vals[cond].append(exp)
                    bias_vals[cond].append(bias)
                    m_vals[cond].append(val)
                    if c == 0:
                        differential_treg_tcon.append(deseq_pval_mat_2[i1, i2] < .2)
                        inds.append(f'{i1}_{i2}')
                else:
                    pass

    inds = arr(inds)
    for cond in exp_vals:
        exp_vals[cond] = arr(exp_vals[cond])
    for cond in bias_vals:
        bias_vals[cond] = arr(bias_vals[cond])
    for cond in m_vals:
        m_vals[cond] = arr(m_vals[cond])
        
    return m_vals, exp_vals, bias_vals, inds, differential_treg_tcon


def process_output(bias_vals, m_vals, exp_vals, inds, filename='distal_benoist_interactions'):
    bias_df = pd.DataFrame(bias_vals)
    bias_df.index = inds
    m_df = pd.DataFrame(m_vals)
    m_df.index = inds
    exp_df = pd.DataFrame(exp_vals)
    exp_df.index = inds

    norm_bias_df = (bias_df.T/bias_df.mean(axis=1)).T
    norm_exp_df = (exp_df.T/exp_df.mean(axis=1)).T

    final_df_benoist = m_df*norm_bias_df
    final_df_benoist = final_df_benoist[~final_df_benoist.isna().any(axis=1)]
    final_df_benoist.astype(int).to_csv(f'./hichip_deseq/{filename}.csv', sep='\t')
    return final_df_benoist

def get_exp_from_indmat(df, chrom, indmat, goodinds):
    exp_mat = np.zeros(indmat.shape)*np.nan
    tmp = df[df['region'] == chrom].set_index('diag')['count.avg']
    ind = np.ravel(indmat[goodinds])
    v = tmp.get(ind).values
    exp_mat[goodinds] = v
    
    
    return exp_mat


def get_bal_exp_from_indmat(df, chrom, indmat, goodinds):
    exp_mat = np.zeros(indmat.shape)*np.nan
    tmp = df[df['region'] == chrom].set_index('diag')['balanced.avg']
    ind = np.ravel(indmat)
    v = tmp.get(ind).values
    exp_mat[ind] = v

    return exp_mat





