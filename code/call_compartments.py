import cooler
import numpy as np
import pandas as pd
import bbi
from copy import deepcopy

arr = np.asarray

import sklearn
from sklearn.decomposition  import PCA

def del_if_chr(ss):
    news = list()
    for s in ss:
        new = s.replace("chr", '')
        news.append(new)
    return news

def compute_chrom_compartment(chrom, cooldict, cond, inds_to_region, chrom_to_start, chrom_to_end, 
                              ignore_sign=False, **kwargs):
    oe_dict, oe_dict_with_pc = make_oe_dicts(cooldict, chrom, condsoi = [cond])
    c, s, e = list(zip(*inds_to_region[chrom]))

    _, pc_treg_250kb, full_pc_treg = make_corr_mat(oe_dict_with_pc[cond], **kwargs)
    if ignore_sign:
        pass
    else:
        h3k27ac_vals = np.ravel(bbi.open('./compartments/treg_h3k27ac.bw').stackup(del_if_chr(c), s, e, bins=1))
        treg_sign = np.sign(np.nanmean(h3k27ac_vals[pc_treg_250kb > 0]) - np.nanmean(h3k27ac_vals[pc_treg_250kb < 0]))
        pc_treg_250kb *= treg_sign

    s = chrom_to_start[chrom]
    e = chrom_to_end[chrom]
    fac = np.sqrt(e-s)
    return pc_treg_250kb*fac, full_pc_treg


def calculate_my_pcs(cooldict, all_ind_to_region, inds_to_region, chrom_to_start, chrom_to_end, 
                     ignore_sign = False,
                     **kwargs):
    compdict = {}
    for cond in cooldict.keys():
        compvec = list()
        seen = set()
        for reg in all_ind_to_region:
            chrom = reg[0]
            if chrom in seen:
                continue
            else:
                seen.add(chrom)
            treg_pc, treg_pca = compute_chrom_compartment(chrom, cooldict, cond, inds_to_region, 
                                                                            chrom_to_start, chrom_to_end, 
                                                                            ignore_sign=ignore_sign, **kwargs)
            compvec += list(treg_pc)
            #full_treg_pcs += list(treg_pca)
            #full_tcon_pcs += list(tcon_pca)
        compdict[cond] = arr(compvec)
    return compdict


def make_corr_mat(m, n_components=2):
    mat = deepcopy(m)
    bad = ((np.isnan(mat)).mean(axis=0) > .5)

    submat = mat[~bad, :][:, ~bad]
    submat[np.isnan(submat)] = 0
    submat[np.isinf(submat)] = 0
    
    pc_input = np.corrcoef(submat)
    pc = PCA(n_components=n_components).fit(pc_input)
#     pc_loadings = PCA(n_components=2).fit_transform(pc_input)
#     pc_loading2 = pc_input@pc.components_[0, :]
    corrs = np.zeros(mat.shape)
    pcs = np.zeros(mat.shape[0])
    pcs[~bad] = pc.components_[0, :]
#     pcs[~bad] = pc_loadings[:, 0]
#     pcs[~bad] = pc_loading2[:, 0]
    pcs[bad] = np.nan
    pc2 = np.zeros(mat.shape[0])
    pc2[~bad] = pc.components_[1, :]
    pc2[bad] = np.nan
    corrs[np.ix_(~bad, ~bad)] = np.corrcoef(submat)  
    return corrs, pcs, (pc.components_, pc2)


def compute_exp(m):
    n = m.shape[1]
    expmat = np.zeros((n, n))
    for i in range(n):
        v = np.nanmean(np.diag(m, k=i))
        np.fill_diagonal(expmat[i:, :], v)
        np.fill_diagonal(expmat[:, i:], v)
    return expmat

def compute_zscore(m):
    n = m.shape[1]
    meanmat = np.zeros((n, n))
    stdmat = np.zeros((n, n))

    for i in range(n):
        mean = np.nanmean(np.diag(m, k=i))
        std = np.nanstd(np.diag(m, k=i))
        np.fill_diagonal(meanmat[i:, :], mean)
        np.fill_diagonal(meanmat[:, i:], mean)
        np.fill_diagonal(stdmat[i:, :], std)
        np.fill_diagonal(stdmat[:, i:], std)

    return meanmat, stdmat


def compute_min_pseudocount(m):
    n = m.shape[1]
    expmat = np.zeros((n, n))
    for i in range(n):
        tmp = np.diag(m, k=i)
        z = tmp[tmp > 0]
        if len(z) == 0:
            v = 0
        else:
            v = np.nanmin(z)
        np.fill_diagonal(expmat[i:, :], v)
        np.fill_diagonal(expmat[:, i:], v)
    return expmat

def make_oe_dicts(cool_dict_big, chrom, plot=True, condsoi=[]):
    oe_dict = {}
    oe_dict_with_pc = {}
    for cond, cool in cool_dict_big.items():
        if len(condsoi) > 0:
            if cond not in condsoi:
                continue
        a = cool.matrix().fetch((chrom))
        a_exp = compute_exp(a)
        pc = compute_min_pseudocount(a)
        oe = np.log2((a)/(a_exp))
        oe_with_pc = np.log2((a+2*pc)/(a_exp+2*pc))

        oe_dict[cond] = oe
        oe_dict_with_pc[cond] = oe_with_pc    
    return oe_dict, oe_dict_with_pc


