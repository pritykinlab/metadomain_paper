from hic_zscore_functions import make_obs_exp_intra, make_obs_exp_inter
from make_figure4 import normalize_raw_intra, normalize_raw_inter
from aux_functions import *

### OLD VERSION
# def construct_oe_mat_OLD(cool, all_ind_to_region, parsed_chroms, chrom_to_start, chrom_to_end, pc=1e-4, log=True):
#     n = len(all_ind_to_region)
#     sep_oe_mat_treg = np.zeros((n, n))
#     blacklist_inds = np.zeros(sep_oe_mat_treg.shape[0])
#     if 'chr' in cool.chromnames[0]:
#         use_chrom = True
#     else:
#         use_chrom = False
#     for i in parsed_chroms:
#         for j in parsed_chroms:
#             s_i, e_i = chrom_to_start[i], chrom_to_end[i]
#             s_j, e_j = chrom_to_start[j], chrom_to_end[j]
#             if use_chrom and 'chr' not in i and 'chr' not in j:
#                 chri = 'chr' + i
#                 chrj = 'chr' + j
#             else:
#                 chri = i
#                 chrj = j
#             if j != i:
#                 m = cool.matrix().fetch(chri, chrj)
#                 f = make_obs_exp_inter(m, nan_to_zero=False, pc=pc, log=log)
#                 sep_oe_mat_treg[s_i:e_i, s_j:e_j] = f
#             elif j == i:
#                 m = cool.matrix().fetch(chri, chrj)
#                 f = make_obs_exp_intra(m, pc=pc, log=log)[0]
#                 if log:
#                     f = np.log2(f)
#                 sep_oe_mat_treg[s_i:e_i, s_j:e_j] = f
#                 blacklist_inds[s_i:e_i] += np.isnan(m).mean(axis=1) > .8
#     blacklist_inds = blacklist_inds > 0
#     return sep_oe_mat_treg, blacklist_inds


def fetch_mat(cool, chrom1, chrom2, padding_size, pc=1e-4, fetch_oe=False, log=True):
    mat = cool.matrix().fetch(chrom1, chrom2)
    if chrom1 == chrom2:
        oe_mat = normalize_raw_intra(mat, nan_to_zero=False, pc = pc, log=log)
    else:
        oe_mat = normalize_raw_inter(mat, nan_to_zero=False, pc = pc, log=log)
    if fetch_oe:
        padded_mat = np.pad(oe_mat, pad_width=padding_size, mode='constant', constant_values=np.nan)
    else:
        padded_mat = np.pad(mat, pad_width=padding_size, mode='constant', constant_values=np.nan)
    return padded_mat

def construct_oe_mat(cool, all_ind_to_region, parsed_chroms, chrom_to_start, chrom_to_end, pc=1e-4, log=True):
    n = len(all_ind_to_region)
    sep_oe_mat_treg = np.zeros((n, n))
    blacklist_inds = np.zeros(sep_oe_mat_treg.shape[0])
    if 'chr' in cool.chromnames[0]:
        use_chrom = True
    else:
        use_chrom = False
    for chrom1 in parsed_chroms:
        for chrom2 in parsed_chroms:
            s_i, e_i = chrom_to_start[chrom1], chrom_to_end[chrom1]
            s_j, e_j = chrom_to_start[chrom2], chrom_to_end[chrom2]
            if use_chrom and 'chr' not in chrom1 and 'chr' not in chrom2:
                chri = 'chr' + chrom1
                chrj = 'chr' + chrom2
            else:
                chri = chrom1
                chrj = chrom2
            m = fetch_mat(cool, chri, chrj, padding_size=0, fetch_oe=True, pc=pc, log=log)
            sep_oe_mat_treg[s_i:e_i, s_j:e_j] = m
            blacklist_inds[s_i:e_i] += np.isnan(m).mean(axis=1) > .8
    blacklist_inds = blacklist_inds > 0
    return sep_oe_mat_treg, blacklist_inds

    

def construct_oe_mat_dict(cooldict, all_ind_to_region, parsed_chroms, chrom_to_start, chrom_to_end, verbose=True, pc=1e-4, log=True):
    oe_mat_dict = {}
    n = len(all_ind_to_region)
    all_blacklist_inds = np.zeros(n)
    for name, cool in cooldict.items():
        oe_mat, blacklist_inds = construct_oe_mat(cool, all_ind_to_region, parsed_chroms, chrom_to_start, chrom_to_end, pc=pc, log=log)
        all_blacklist_inds += blacklist_inds
        oe_mat_dict[name] = oe_mat
        if verbose:
            print("Done with", name)
    all_blacklist_inds = all_blacklist_inds > 0
    for name, oe_mat in oe_mat_dict.items():
        oe_mat[all_blacklist_inds, :] = np.nan
        oe_mat[:, all_blacklist_inds] = np.nan
    return oe_mat_dict