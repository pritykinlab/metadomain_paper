import scipy
import scipy.signal
import numpy as np
import pandas as pd
from hic_zscore_functions import make_obs_exp_nolog
from copy import deepcopy
arr = np.asarray

def make_auxiliary_scores(chrom1, chrom2, mat, oe_mat_ind, A_mat_ind, B_mat_ind, chrom_to_start_50kb):    
    n1, n2 = mat.shape
    inner_scores = []
    outer_scores = []
    vert_scores = []
    side_scores = []

    c = 0
    inner_filter, outer_filter, left_filter, right_filter, up_filter, down_filter = set_filters()
    s = inner_filter.shape[0]//2
    side_filter = left_filter + right_filter
    vert_filter = up_filter + down_filter

    offset = s
    filters = [outer_filter, outer_filter, outer_filter]
    matrices = [A_mat_ind*mat, B_mat_ind*mat, oe_mat_ind*mat]
    labels = ['outer_A_scores', 'outer_B_scores', 'outer_oe_scores']
    score_df = process_matrix_with_filters(chrom1, chrom2, matrices, filters, labels, chrom_to_start_50kb)
    return score_df


def set_filters():
    n = 31
    inner_filter = np.zeros((n, n))
    outer_filter = np.ones((n, n))

    w = 3
    mid = inner_filter.shape[1]//2
    inner_filter[mid-w:mid+w+1, mid-w:mid+w+1] = 1
    outer_filter -= inner_filter

    left_filter = np.zeros((n, n))
    left_filter[mid-w:mid+w+1, :n//2] = 1
    left_filter[inner_filter==1] = 0

    right_filter = np.zeros((n, n))
    right_filter[mid-w:mid+w+1, n//2:] = 1
    right_filter[inner_filter==1] = 0

    up_filter = left_filter.T
    down_filter = right_filter.T
    return inner_filter, outer_filter, left_filter, right_filter, up_filter, down_filter

def make_ind_dict(chrom, merged_cool, compartment_vec, chrom_to_start, chrom_to_end):
    inddict = {}
    m = merged_cool.matrix(balance=True).fetch(chrom).astype(float)
    oe_mat = np.log2(make_obs_exp_nolog(m, pc=0)[0])
    inddict["oe"] = oe_mat > 0
    s, e = chrom_to_start[chrom], chrom_to_end[chrom]
    comp = deepcopy(compartment_vec)[s:e]
    assert len(comp) == len(m)
    A_comp_ind = np.add.outer(comp, comp) > 1.2
    B_comp_ind = np.add.outer(comp, comp) < -1.2
    inddict["A"] = A_comp_ind
    inddict["B"] = B_comp_ind
    return inddict

def make_count_df(chrom1, chrom2, inddict, chrom_to_start):
    countdict = {}
    inner_filter, outer_filter, left_filter, right_filter, up_filter, down_filter = set_filters()
    oe_mat_ind, A_mat_ind, B_mat_ind = inddict['oe'], inddict['A'], inddict['B']
    filters = [outer_filter, outer_filter, outer_filter]
    matrices = [A_mat_ind, B_mat_ind, oe_mat_ind]
    labels = ['outer_A_counts', 'outer_B_counts', 'outer_oe_counts']
    count_df = process_matrix_with_filters(chrom1, chrom2, matrices, filters, labels, chrom_to_start)
    return count_df 
    
def prepare_raw_matrix(cool, chrom):
    m_tmp = np.triu(cool.matrix(balance=False).fetch(chrom).astype(float))
    np.fill_diagonal(m_tmp, np.nan)
    np.fill_diagonal(m_tmp[1:, :], np.nan)
    np.fill_diagonal(m_tmp[:, 1:], np.nan)
    return m_tmp

def prepare_raw_inter_matrix(cool, chrom1, chrom2):
    m_tmp = cool.matrix(balance=False).fetch(chrom1, chrom2).astype(float)
    return m_tmp


def get_aux_df(chrom1, chrom2, cond, m_tmp, cool, inddict, countdict, chrom_to_start,):    
    oe_mat_ind, A_comp_ind, B_comp_ind = inddict["oe"], inddict["A"], inddict["B"]
    aux_score_df = make_auxiliary_scores(chrom1, chrom2, m_tmp, oe_mat_ind, A_comp_ind, B_comp_ind, chrom_to_start)
    return aux_score_df


def get_oe_mat(mat, pc=0):
    oe_mat = np.log2(make_obs_exp_nolog(mat, pc=0)[0])
    return oe_mat

def get_inter_exp_mat(mat, pc=0):
    e = np.nanmean(mat)
    return e

def get_inter_oe_mat(mat, pc=0):
    e = np.nanmean(mat)
    oe_mat = np.log2((mat+pc)/(e+pc))
    return oe_mat


def get_places_with_high_oe(oe_mat, chrom1, chrom2, chrom_to_start,  co=0):
    indsoi = arr(np.where(oe_mat > co))
    X, Y = indsoi 
    X = X + chrom_to_start[chrom1]
    Y = Y + chrom_to_start[chrom2]
    v = np.ravel(X).astype(str)
    w = np.ravel(Y).astype(str)
    places_with_high_oe = pd.Series(v) + "_" + pd.Series(w)
    return places_with_high_oe

def normalize_aux_df(aux_df, count_df, cond, tmp_df):
    tmp_df["A_scores_" + cond] = ((aux_df['outer_A_scores'] / count_df['outer_A_counts'])*10).astype(float).round()
    tmp_df["B_scores_" + cond] = ((aux_df['outer_B_scores'] / count_df['outer_B_counts'])*10).astype(float).round()
    tmp_df["oe_scores_" + cond] = ((aux_df['outer_oe_scores'] / count_df['outer_oe_counts'])*10).astype(float).round()
    tmp_df["A_scores_" + cond] = tmp_df["A_scores_" + cond].fillna(0)
    tmp_df["B_scores_" + cond] = tmp_df["B_scores_" + cond].fillna(0)
    tmp_df["oe_scores_" + cond] = tmp_df["oe_scores_" + cond].fillna(0)    
    return tmp_df



def make_readcount_df(params, n_filters = 4):
    (merged_cool, chrom_to_start, 
      chrom_to_end, cool_dict, chrom, compartment_vec) = params

    balanced_mat = merged_cool.matrix(balance=True).fetch(chrom).astype(float)
    nanmarker = np.isnan(balanced_mat)
    oe_mat = np.triu(get_oe_mat(balanced_mat), k = 32)
    places_with_high_oe = get_places_with_high_oe(oe_mat, chrom, chrom, chrom_to_start, co=1)
    nan_df = make_scores(chrom, chrom, np.triu(~nanmarker), chrom_to_start)

    t = nan_df.index.isin(places_with_high_oe)
    nan_df = nan_df.loc[t]
    inddict = make_ind_dict(chrom, merged_cool, compartment_vec,
                           chrom_to_start, chrom_to_end)
    count_df = make_count_df(chrom, chrom, inddict, chrom_to_start)
    assert len(nan_df.columns) == 4 ## One column is for 'places'
    frac_nan_df = (nan_df)/(nan_df.max(axis=0))

    tmp_df = pd.DataFrame()
    for cond, cool in cool_dict.items():
        m_tmp = prepare_raw_matrix(cool, chrom)
        score_df = make_scores(chrom, chrom, m_tmp, chrom_to_start)
        score_df = score_df.loc[t]
        score_df = (score_df/frac_nan_df).astype(float).round()
        tmp_df["inner_" + cond] = score_df['inner_scores']
        tmp_df["outer_" + cond] = score_df['outer_scores']
        tmp_df["vert_" + cond] = score_df['vert_scores']
        tmp_df["side_" + cond] = score_df['side_scores']
        #break
    for cond, cool in cool_dict.items():
        m_tmp = prepare_raw_matrix(cool, chrom)
        aux_df = get_aux_df(chrom, chrom, cond, m_tmp, cool, inddict, count_df, chrom_to_start)
        aux_df = aux_df.loc[t]
        sub_count_df = count_df.loc[t]
        normalize_aux_df(aux_df, sub_count_df, cond, tmp_df)
        #break
    tmp_df.index = aux_df.index
    return tmp_df


import scipy
n_workers = 16
def make_scores(chrom1, chrom2, mat, chrom_to_start_50kb):
    n1, n2 = mat.shape

    c = 0
    inner_filter, outer_filter, left_filter, right_filter, up_filter, down_filter = set_filters()
    s = inner_filter.shape[0]//2
    side_filter = left_filter + right_filter
    vert_filter = up_filter + down_filter

    offset = s
    filters = [inner_filter, outer_filter, vert_filter, side_filter]
    labels = ['inner_scores', 'outer_scores', 'vert_scores', 'side_scores']
    matrices = (mat, mat, mat, mat)
    score_df = process_matrix_with_filters(chrom1, chrom2, matrices, filters, labels, chrom_to_start_50kb)
    return score_df

import scipy
import scipy.signal
def process_matrix_with_filters(chrom1, chrom2, matrices, filters, labels, chrom_to_start_50kb):    
    score_df = pd.DataFrame()
    for filt, label, mat in zip(filters, labels, matrices):
        n1, n2 = mat.shape
        c = 0
        s = filt.shape[0]//2
        offset = s
        scores = scipy.signal.correlate2d(mat, filt)[offset:-offset, offset:-offset]
        score_df[label] = np.ravel(scores)

    X, Y = np.indices(mat.shape)
    X = X + chrom_to_start_50kb[chrom1]
    Y = Y + chrom_to_start_50kb[chrom2]
    v = np.ravel(X).astype(str)
    w = np.ravel(Y).astype(str)
    places = pd.Series(v) + "_" + pd.Series(w)        
    score_df.index = np.ravel(places)
    return score_df


def make_interchrom_ind_dict(chrom1, chrom2, merged_cool, compartment_vec, chrom_to_start, chrom_to_end):
    inddict = {}
    m = merged_cool.matrix(balance=True).fetch(chrom1, chrom2).astype(float)
    oe_mat = get_inter_oe_mat(m)
    inddict["oe"] = oe_mat > 0
    s1, e1 = chrom_to_start[chrom1], chrom_to_end[chrom1]
    s2, e2 = chrom_to_start[chrom2], chrom_to_end[chrom2]
    comp1 = deepcopy(compartment_vec)[s1:e1]
    comp2 = deepcopy(compartment_vec)[s2:e2]
    assert len(comp1) == m.shape[0]
    assert len(comp2) == m.shape[1]
    A_comp_ind = np.add.outer(comp1, comp2) > 1.2
    B_comp_ind = np.add.outer(comp1, comp2) < -1.2
    inddict["A"] = A_comp_ind
    inddict["B"] = B_comp_ind
    return inddict


def make_interchrom_readcount_df(params, n_filters = 4):
    (merged_cool, chrom_to_start,
      chrom_to_end, cool_dict, chrom1, chrom2, compartment_vec) = params

    balanced_mat = merged_cool.matrix(balance=True).fetch(chrom1, chrom2).astype(float)
    nanmarker = np.isnan(balanced_mat)
    oe_mat = get_inter_oe_mat(balanced_mat)
    places_with_high_oe = get_places_with_high_oe(oe_mat, chrom1, chrom2, chrom_to_start, co=0)

    nan_df = make_scores(chrom1, chrom2, ~nanmarker, chrom_to_start)
    t = nan_df.index.isin(places_with_high_oe)
    nan_df = nan_df.loc[t]
    inddict = make_interchrom_ind_dict(chrom1, chrom2, merged_cool, compartment_vec,
                           chrom_to_start, chrom_to_end)
    count_df = make_count_df(chrom1, chrom2, inddict, chrom_to_start)
    assert len(nan_df.columns) == 4 ## One column is for 'places'
    frac_nan_df = (nan_df)/(nan_df.max(axis=0))

    tmp_df = pd.DataFrame()
    for cond, cool in cool_dict.items():
        m_tmp = prepare_raw_inter_matrix(cool, chrom1, chrom2)
        score_df = make_scores(chrom1, chrom2, m_tmp, chrom_to_start)
        score_df = score_df.loc[t]
        score_df = (score_df/frac_nan_df).astype(float).round()
        tmp_df["inner_" + cond] = score_df['inner_scores']
        tmp_df["outer_" + cond] = score_df['outer_scores']
        tmp_df["vert_" + cond] = score_df['vert_scores']
        tmp_df["side_" + cond] = score_df['side_scores']
        #break
    for cond, cool in cool_dict.items():
        m_tmp = prepare_raw_inter_matrix(cool, chrom1, chrom2)
        aux_df = get_aux_df(chrom1, chrom2, cond, m_tmp, cool, inddict, count_df, chrom_to_start)
        aux_df = aux_df.loc[t]
        sub_count_df = count_df.loc[t]
        normalize_aux_df(aux_df, sub_count_df, cond, tmp_df)
        #break
    tmp_df.index = aux_df.index
    return tmp_df, chrom1, chrom2, t









