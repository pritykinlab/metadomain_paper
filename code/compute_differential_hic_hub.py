import statsmodels
import pandas as pd
import statsmodels.stats
from aux_functions import *
import statsmodels.stats.multitest
from hub_pileup import append_mats

def make_outside_filt(m, center, method = 'center_only', square=None):
    if square is None:
        square = center
    else:
        square = square
    outside_filt = np.ones((2*m+1, 2*m+1))
    outside_filt[m-center : m+center+1, m-center: m+center+1] = 0

    inside_filt = np.zeros((2*m+1, 2*m+1))
    inside_filt[m-center : m+center+1, m-center: m+center+1] = 1

    if (method == 'center_only'):
        outside_filt[:m] = 0
        outside_filt[m+1:] = 0
        inside_filt[:m] = 0
        inside_filt[m+1:] = 0
    elif method == 'center_square':
        outside_filt[:m-square ] = 0
        outside_filt[m+square+1:] = 0
        inside_filt[:m-center] = 0
        inside_filt[m+center+1:] = 0
    elif method == 'center_square_and_skip':
        outside_filt[:m-center ] = 0
        outside_filt[m+center+1:] = 0
        outside_filt[:, m-center-center: m+center+center+1] = 0
        inside_filt[:m-center] = 0
        inside_filt[m+center+1:] = 0

    elif method == 'outside_square_and_skip':
        outside_filt[ m-center-center: m+center+center+1, m-center-center: m+center+center+1] = 0
        # inside_filt[:m-center] = 0
        # inside_filt[m+center+1:] = 0


    outside_filt = outside_filt>0
    inside_filt = inside_filt>0

    return inside_filt, outside_filt


# def make_sideways_filt(m, center, delta):
#     outside_filt = np.zeros((2*m+1, 2*m+1))
#     outside_filt[m-center : m+center+1, :] = 1

#     inside_filt = np.zeros((2*m+1, 2*m+1))
#     inside_filt[m-center : m+center+1, m-center: m+center+1] = 1
#     outside_filt[m-delta : m+delta+1, m-delta: m+delta+1] = 0
    
#     outside_filt[inside_filt>0] = 0
#     outside_filt = outside_filt>0
#     inside_filt = inside_filt>0

#     return inside_filt, outside_filt


def compile_pileup_mat_for_ind_over_hub(ind, hub_inds, sep_oe_mat_treg, outside_filt):
    ms_treg = []
    m = len(outside_filt)//2
    for i in hub_inds:
        ms_treg.append(sep_oe_mat_treg[ind-m:ind+1+m, i-m:i+m+1].copy())
    ms_treg = arr(ms_treg)
    return ms_treg


def get_inside_outside(ms_treg, inside_filt, outside_filt):
    m = len(outside_filt)//2
    ms_treg = arr(ms_treg)
    v_middle = np.nanmean(ms_treg[:, inside_filt], axis=1)
    v_outside = np.nanmean(ms_treg[:, outside_filt], axis=1)
    return v_middle, v_outside


def get_differential_hic(sep_oe_mat_treg, sep_oe_mat_tcon, hub_dict, all_ind_to_region, m = 2, center = 1, cliplo=-1, cliphigh=10,
                         method = 'center_square', meta_inds = {}):

    inside_filt, outside_filt = make_outside_filt(m, center, method=method)
    # print(np.sum(inside_filt), np.sum(outside_filt))
    # plt.matshow(inside_filt)
    # plt.matshow(outside_filt)
    # return inside_filt, outside_filt
    tstat_df = pd.DataFrame()
    stat_df = pd.DataFrame()
    pval_df = pd.DataFrame()
    raw_pval_df = pd.DataFrame()
    meta = [] 
    for u, idx in hub_dict.items():
        xs = np.arange(10, len(all_ind_to_region)-10)
        stats, ps, deltas = [], [], []
        for x in xs:
            if (len(meta_inds) > 0) and (x not in meta_inds):
                continue
            treg_pileup_mat = compile_pileup_mat_for_ind_over_hub(x, idx, sep_oe_mat_treg, outside_filt)
            tcon_pileup_mat = compile_pileup_mat_for_ind_over_hub(x, idx, sep_oe_mat_tcon, outside_filt)
            inside_treg, outside_treg = get_inside_outside(treg_pileup_mat, inside_filt, outside_filt)
            inside_tcon, outside_tcon = get_inside_outside(tcon_pileup_mat, inside_filt, outside_filt)
            stat, p, delta, _, _ = test_inside_outside_differential(inside_treg, inside_tcon, outside_treg, outside_tcon, cliplo=cliplo, cliphigh=cliphigh)
            stats.append(stat)
            ps.append(p)
            deltas.append(delta)
            if (x in meta_inds):
                meta.append([u, x, inside_treg, inside_tcon, outside_treg, outside_tcon])
        new_pval = arr(ps).copy()
        idx = ~np.isnan(new_pval)
        new_pval[idx] = statsmodels.stats.multitest.fdrcorrection(new_pval[idx])[1]
        tstat_df[u] = stats
        stat_df[u] = deltas
        pval_df[u] = new_pval
        raw_pval_df[u] = ps
        if len(meta_inds) > 0:
            pass
        else:
            tstat_df.index = xs
            stat_df.index = xs
            pval_df.index = xs
            raw_pval_df.index = xs
    return tstat_df, stat_df, pval_df, raw_pval_df, meta

def test_inside_outside_differential(v_inside1, v_inside2, v_outside1, v_outside2, cliplo=-1, cliphigh=10, nan_policy='omit',test=scipy.stats.ranksums):
    v_middle = v_inside1 - v_inside2
    v_outside = v_outside1 - v_outside2
    stat, p = test(v_middle, v_outside, nan_policy=nan_policy)
    delta = np.nanmean(v_middle.clip(cliplo, cliphigh)) - np.nanmean(v_outside.clip(cliplo, cliphigh))
    return stat, p, delta, v_middle, v_outside





def test_inside_outside_baseline(v_middle, v_outside, cliplo=-1, cliphigh=10):
    stat, p = scipy.stats.ranksums(v_middle, v_outside, nan_policy='omit')
    delta = np.nanmean(v_middle.clip(cliplo, cliphigh)) - np.nanmean(v_outside.clip(cliplo, cliphigh))
    return stat, p, delta, v_middle, v_outside



import itertools
from make_figure4 import *



# def get_differential_hic_pileup(ind_df, inter_and_intra_connections_treg, 
#                                 inter_and_intra_connections_tcon, chrom_to_start,
#                                 megaloop_pileup_cooldict = {}, padding_size = 30):
#     all_mat_dict = defaultdict(list)
#     all_metadata = []

#     chromsoi = ind_df['chrom'].unique()

#     for chrom1, chrom2 in itertools.product(chromsoi, chromsoi):
#         if chrom1 <= chrom2:
#             continue
#         else:
#             idx1 = ind_df['chrom'] == chrom1
#             idx2 = ind_df['chrom'] == chrom2
#             subdf_X = ind_df[idx1].copy()
#             subdf_Y = ind_df[idx2].copy()
#             subdf_X['norm_ind'] = subdf_X['ind'] - chrom_to_start[chrom1]
#             subdf_Y['norm_ind'] = subdf_Y['ind'] - chrom_to_start[chrom2]
#             rows, cols = list(zip(*itertools.product(subdf_X.index, subdf_Y.index)))
            
#             for name, cool in megaloop_pileup_cooldict.items():
#                 append_mats(subdf_X, subdf_Y, cool, name, rows, cols, all_mat_dict, padding_size=padding_size)

#             for row, col in zip(rows, cols):
#                 row_X = subdf_X.loc[row]
#                 row_Y = subdf_Y.loc[col]

#                 row, col = 5*row_X['norm_ind']+padding_size, 5*row_Y['norm_ind']+padding_size
#                 og_row, og_col = row_X['norm_ind'], row_Y['norm_ind']
#                 ind_X, ind_Y = row_X['ind'], row_Y['ind']
#                 cluster_X, cluster_Y = row_X['cluster'], row_Y['cluster']

#                 mega_loops_treg = inter_and_intra_connections_treg[ind_X, ind_Y]
#                 mega_loops_tcon = inter_and_intra_connections_tcon[ind_X, ind_Y]
#                 all_metadata.append([mega_loops_treg, mega_loops_tcon, ind_X, ind_Y, cluster_X, cluster_Y])

#     metadata_df = pd.DataFrame(all_metadata, columns = ['treg_mega', 'tcon_mega', 'ind1', 'ind2', 'cluster1', 'cluster2'])

#     for key, v in all_mat_dict.items():
#         all_mat_dict[key] = np.stack(all_mat_dict[key])    
#     return all_mat_dict, metadata_df



# def get_ind_differential_hic_pileup(ind_df1, ind_df2, inter_and_intra_connections_treg, 
#                                 inter_and_intra_connections_tcon, chrom_to_start,
#                                 megaloop_pileup_cooldict = {}, padding_size = 30, pc=1e-4, log=True):
#     all_mat_dict = defaultdict(list)
#     all_metadata = []

#     chromsoi = sorted(list(set(np.unique(list(ind_df1['chrom'].unique()) + list(ind_df2['chrom'].unique())))))
#     for chrom1, chrom2 in itertools.product(chromsoi, chromsoi):
#         if chrom1 == chrom2:
#             continue
#         else:
#             idx1 = ind_df1['chrom'] == chrom1
#             idx2 = ind_df2['chrom'] == chrom2
#             subdf_X = ind_df1[idx1].copy()
#             subdf_Y = ind_df2[idx2].copy()
#             if len(subdf_X) == 0 or len(subdf_Y) == 0:
#                 continue

#             subdf_X['norm_ind'] = subdf_X['ind'] - chrom_to_start[chrom1]
#             subdf_Y['norm_ind'] = subdf_Y['ind'] - chrom_to_start[chrom2]
#             rows, cols = list(zip(*itertools.product(subdf_X.index, subdf_Y.index)))
            
#             for name, cool in megaloop_pileup_cooldict.items():
#                 append_mats(subdf_X, subdf_Y, cool, name, rows, cols, all_mat_dict, pc=pc, padding_size=padding_size, log=log)

#             for row, col in zip(rows, cols):
#                 row_X = subdf_X.loc[row]
#                 row_Y = subdf_Y.loc[col]

#                 row, col = 5*row_X['norm_ind']+padding_size, 5*row_Y['norm_ind']+padding_size
#                 ind_X, ind_Y = row_X['ind'], row_Y['ind']
#                 cluster_X, cluster_Y = row_X['cluster'], row_Y['cluster']

#                 mega_loops_treg = inter_and_intra_connections_treg[ind_X, ind_Y]
#                 mega_loops_tcon = inter_and_intra_connections_tcon[ind_X, ind_Y]
#                 all_metadata.append([mega_loops_treg, mega_loops_tcon, ind_X, ind_Y, cluster_X, cluster_Y])

#     metadata_df = pd.DataFrame(all_metadata, columns = ['treg_mega', 'tcon_mega', 'ind1', 'ind2', 'cluster1', 'cluster2'])

#     for key, v in all_mat_dict.items():
#         all_mat_dict[key] = np.stack(all_mat_dict[key])    
#     return all_mat_dict, metadata_df



# def get_inside(ind, hub_inds, sep_oe_mat_treg, outside_filt):
#     ms_treg = []
#     ms_tcon = []
#     m = len(outside_filt)//2
#     for i in hub_inds:
#         ms_treg.append(sep_oe_mat_treg[ind-m:ind+1+m, i-m:i+m+1].copy())
#     ms_treg = arr(ms_treg)
#     ms_tcon = arr(ms_tcon)
#     if ms_treg.shape[1:] != (2*m+1, 2*m+1):
#         stat = np.nan
#         p = np.nan
#         delta = np.nan  
#         v_middle = np.nan
#         v_outside = np.nan
#     else:
#         v_middle1 = ms_treg[:, m, m]
#     return v_middle1

# def get_outside(ind, hub_inds, sep_oe_mat_treg, outside_filt):
#     ms_treg = []
#     ms_tcon = []
#     m = len(outside_filt)//2
#     for i in hub_inds:
#         ms_treg.append(sep_oe_mat_treg[ind-m:ind+1+m, i-m:i+m+1].copy())
#     ms_treg = arr(ms_treg)
#     ms_tcon = arr(ms_tcon)
#     if ms_treg.shape[1:] != (2*m+1, 2*m+1):
#         stat = np.nan
#         p = np.nan
#         delta = np.nan  
#         v_middle = np.nan
#         v_outside = np.nan
#     else:
#         v_middle1 = np.nanmean(ms_treg[:, outside_filt], axis=1)
#     return v_middle1


# def fetch_mat(chrom1, chrom2, cool, fetch_raw=False, padding_size = 30, pc = 1e-4, log=True):
#     if 'chr' not in chrom1:
#         mat = cool.matrix().fetch('chr' + chrom1, 'chr' + chrom2)
#     else:
#         mat = cool.matrix().fetch(chrom1, chrom2)
#     if chrom1 == chrom2:
#         oe_mat = make_obs_exp_intra(mat, pc=pc, log=log)[0]
#     else:
#         oe_mat = make_obs_exp_inter(mat, nan_to_zero=False, pc=pc, log=log)
#     if fetch_raw:
#         padded_mat = np.pad(mat, pad_width=padding_size, mode='constant', constant_values=np.nan)
#     else:
#         padded_mat = np.pad(oe_mat, pad_width=padding_size, mode='constant', constant_values=np.nan)
#     return padded_mat
    

# def append_mats(subdf_X, subdf_Y, cool, name, rows, cols, mat_dict, fetch_raw=False, padding_size = 30, pc=1e-4, stride=5,
#                 sl1=16, sl2=-13):
#     chrom1, chrom2 = subdf_X['chrom'].iloc[0], subdf_Y['chrom'].iloc[0]
#     padded_mat = fetch_mat(chrom1, chrom2, cool, fetch_raw=fetch_raw, padding_size = padding_size, pc=pc)
#     for row, col in zip(rows, cols):
#         row_X = subdf_X.at[row, 'norm_ind']
#         row_Y = subdf_Y.at[col, 'norm_ind']
#         row, col = stride*row_X+padding_size, stride*row_Y+padding_size
#         submat = padded_mat[row-padding_size:row+padding_size, 
#                             col-padding_size:col+padding_size][sl1:sl2, sl1:sl2]
#         # print(submat.shape, padded_mat.shape, row, col)
#         # assert submat.shape == (31, 31), (submat.shape, row, col, padded_mat.shape)
#         mat_dict[name].append(submat)




# def test_inside_outside_differential(ind, hub_inds, sep_oe_mat_treg, sep_oe_mat_tcon, outside_filt, cliplo=-1, cliphigh=10):
#     ms_treg = []
#     ms_tcon = []
#     m = len(outside_filt)//2
#     for i in hub_inds:
#         ms_treg.append(sep_oe_mat_treg[ind-m:ind+1+m, i-m:i+m+1].copy())
#         ms_tcon.append(sep_oe_mat_tcon[ind-m:ind+1+m, i-m:i+m+1].copy())
#     ms_treg = arr(ms_treg)
#     ms_tcon = arr(ms_tcon)
#     if ms_treg.shape[1:] != (2*m+1, 2*m+1):
#         stat = np.nan
#         p = np.nan
#         delta = np.nan  
#         v_middle = np.nan
#         v_outside = np.nan
#     else:
#         v_middle = ms_treg[:, m, m] - ms_tcon[:, m, m]
#         v_outside = np.nanmean(ms_treg[:, outside_filt], axis=1) - np.nanmean(ms_tcon[:, outside_filt], axis=1)
#         stat, p = scipy.stats.ranksums(v_middle, v_outside, nan_policy='omit')
#         delta = np.nanmean(v_middle.clip(cliplo, cliphigh)) - np.nanmean(v_outside.clip(cliplo, cliphigh))
#     return stat, p, delta, v_middle, v_outside


