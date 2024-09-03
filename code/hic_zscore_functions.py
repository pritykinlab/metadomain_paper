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
from copy import deepcopy


def make_expected_for_zscore(raw_mat):
    means = []
    for i in range(len(raw_mat)):
        m = np.nanmean(np.diag(raw_mat, k=i))
        if np.isnan(m):
            break
        means.append(m)
    print("Max diag for expected zscore", i)
    meanmat = np.zeros(raw_mat.shape)
    diags = np.triu(np.flip(np.indices(raw_mat.shape)[0]) + (np.indices(raw_mat.shape)[1]))
    diags += (diags).T - np.diag(np.diag(diags))-diags[0, 0]
    exp = deepcopy(diags).astype(float)
    exp[np.isnan(raw_mat)] = np.nan
    for c, val in enumerate(means):
        np.fill_diagonal(exp[c:, :], val)
        np.fill_diagonal(exp[:, c:], val)
        
    stds = []
    for i in range(len(raw_mat)):
        m = np.nanstd(np.diag(raw_mat, k=i))
        if np.isnan(m):
            stds.append(np.nan)
            break
        elif i == 0:
            stds.append(m)
        elif m == 0:
            stds.append(stds[-1])
        else:
            if ~np.isnan(m):
                stds.append(m)
            else:
                stds.append(stds[-1])
    diags = np.triu(np.flip(np.indices(raw_mat.shape)[0]) + (np.indices(raw_mat.shape)[1]))
    diags += (diags).T - np.diag(np.diag(diags))-diags[0, 0]
    std_mat = deepcopy(diags).astype(float)
    std_mat[np.isnan(std_mat)] = np.nan
    for c, val in enumerate(stds):
        np.fill_diagonal(std_mat[c:, :], val)
        np.fill_diagonal(std_mat[:, c:], val)
        if np.isnan(val):
            std_mat[diags >= c] = val
    return exp, std_mat

def make_diag_zscore_no_mean(mat):
    mat[np.isinf(mat)] = 0 
    exp, std = make_expected_for_zscore(mat)
    return (mat - 0)/std, exp, std

def make_diag_zscore_smooth(mat):
    mat[np.isinf(mat)] = 0 
    exp, std = make_expected_for_zscore(mat)
    xs = np.log10(np.arange(1, 1+len(std[0, :])))
    ys = np.log10(deepcopy(std[0, :]))

    bad = np.where(np.isnan(ys))[0][0]
    ys = ys[:bad]
    xs = xs[:bad]

    fit = np.poly1d(np.polyfit(xs, ys, 3))
    smoothed_stds = np.power(10, fit(xs))
    # last_ind = np.where(~np.isnan(smoothed_stds))[0][-1]
    # print(last_ind, smoothed_stds[last_ind], smoothed_stds[last_ind-1])
    # smoothed_stds[last_ind] = smoothed_stds[last_ind-1]

    for c, val in enumerate(smoothed_stds):
        np.fill_diagonal(std[c:, :], val)
        np.fill_diagonal(std[:, c:], val)
    return (mat - exp)/std, exp, std


def make_diag_zscore(mat):
    mat[np.isinf(mat)] = 0 
    exp, std = make_expected_for_zscore(mat)
    return (mat - exp)/std, exp, std

def make_trans_zscore(mat, std=None):
    mat[np.isinf(mat)] = 0 
    if std:
        exp, _ = np.nanmean(mat), np.nanstd(mat) 
    else:
        exp, std = np.nanmean(mat), np.nanstd(mat)
    return (mat - exp)/std, exp, std

def make_obs_exp_nolog(balanced_mat, mat_type='intra', pc=1e-6):
    exp = make_expected(balanced_mat, mat_type=mat_type)
    f = (balanced_mat+pc)/(exp+pc)
    return f, exp

# def make_obs_exp_balanced(raw_mat, pc=1e-6, balance=False, nan_to_zero=True):
#     exp = make_expected(raw_mat)
#     f = np.log2((raw_mat+pc)/(exp+pc))
#     f[np.isinf(f)] = 0
#     if nan_to_zero:
#         f[np.isnan(f)] = 0
#     return f, exp

def make_obs_exp_intra(balanced_mat, pc=1e-6, log=False, nan_to_zero=False):
    assert is_symmetric(balanced_mat)
    exp = make_expected(balanced_mat, mat_type='intra')
    f = (balanced_mat+pc)/(exp+pc)
    if log:
        f = np.log2(f)
    if nan_to_zero:
        f[np.isnan(f)] = 0
    return f, exp


def make_obs_exp_inter(raw_mat, pc=1e-8,  nan_to_zero=False, log=True):
    assert not is_symmetric(raw_mat)
    exp = make_expected(raw_mat, mat_type='inter')
    f = (raw_mat+pc)/(exp+pc)
    if log:
        f = np.log2(f)
    if nan_to_zero:
        f[np.isnan(f)] = 0
    return f    





def is_symmetric(mat, tol=1e-8):
    if mat.shape[0] != mat.shape[1]:
        return False
    return np.nansum(np.abs(mat-mat.T)) < tol
    
def make_expected(balanced_mat, mat_type='intra'):
    if (mat_type == 'inter') or (not is_symmetric(balanced_mat)):
        return np.nanmean(balanced_mat)
    
    exp = np.zeros(balanced_mat.shape)
    iu_x, iu_y = np.diag_indices(len(balanced_mat))
    last_m = 0
    for off_diag_k in range(len(balanced_mat)):
        m = np.nanmean(np.diag(balanced_mat, k=off_diag_k))
        if off_diag_k > 0 and np.isnan(m):       
            exp[(iu_x, iu_y)] = last_m 
            exp[(iu_y, iu_x)] = last_m 
        else:
            exp[(iu_x, iu_y)] = m
            exp[(iu_y, iu_x)] = m
        iu_x, iu_y = iu_x[:-1], iu_y[1:]
        last_m = m
    return exp


