import numpy as np
import scipy


# def filt(x):
#     return scipy.ndimage.gaussian_filter1d(x, sigma=12)

def make_contained(hicdc_ancs, hicdc_loops, mustache_full_anchors, mustache_full_loops, slopsize=30000):
    slop_to_unslopped = {}
    for i, j in zip(hicdc_ancs.slop(b=slopsize, g='./annotations/chromsizes'), hicdc_ancs):
        l1, l2 = tuple(i[:3]), tuple(j[:3])
        slop_to_unslopped[l1] = l2


    hicdc_loops_dict = {}
    for i in hicdc_loops:
        l1, l2 = tuple(i[:3]), tuple(i[3:6])
        hicdc_loops_dict[(l1, l2)]= 1
        
    mustache_to_hicdc = {}
    for i in mustache_full_anchors:
        mustache_anc = tuple(i[:3])
        mustache_to_hicdc[mustache_anc] = []

    for i in hicdc_ancs.slop(b=slopsize, g='./annotations/chromsizes').intersect(mustache_full_anchors, wo=True):
        mustache_anc = tuple(i[3:6])
        hicdc_anc = slop_to_unslopped[tuple(i[:3])]
        mustache_to_hicdc[mustache_anc].append(hicdc_anc)

        
    contained = []
    not_contained = []
    c = 0
    for i in mustache_full_loops:
        temp=0
        l1, l2 = tuple(i[:3]), tuple(i[3:6])
        for i in mustache_to_hicdc[l1]:
            for j in mustache_to_hicdc[l2]:    
                if hicdc_loops_dict.get((i, j)):
                    temp  = 1
        if temp == 0:
            not_contained.append(list(l1)+list(l2))
        if temp == 1:
            contained.append(list(l1)+list(l2))
        c += temp
    print(c/len(mustache_full_loops))
    return contained, not_contained                


def do_interp(y):
    nans, x= nan_helper(y)
    if np.nan_to_num(y)[10] < np.nan_to_num(y)[-10]:
        y[nans]= np.interp(x(nans), x(~nans), y[~nans], left=0)
    else:
        y[nans]= np.interp(x(nans), x(~nans), y[~nans], right=0)        
    return y


def merge_mask(mask):
    mask_where = np.where(mask> 0)[0]
    newmask = deepcopy(mask)
    mergesize = 20
    for i, _ in enumerate(mask_where[1:]):
        if (mask_where[i] -mask_where[i-1] > 1) and (mask_where[i] - mask_where[i-1] < mergesize):
            newmask[mask_where[i-1]:mask_where[i+1]] = 1
    return newmask


def sum_over_mask(mask, r):
    r = np.expand_dims(r, axis=0)
    labels, unique_count = scipy.ndimage.label(mask)
    unique_vals = np.unique(labels)[1:]
    # mat = np.zeros((len(mask), unique_count))
    newr = np.zeros(r.shape)
    for _, u in enumerate(unique_vals):
        newvec = np.zeros(len(mask))
        newvec[labels==u] = 1/(labels==u).sum()
        newr[:, newvec >0] = np.expand_dims(r@newvec, axis=1)    
    return newr.reshape(newr.shape[1])

def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]


def hicdc_do_calculation(chrom, s, e, treg_expecteds_balanced, tcon_expecteds_balanced, cool_treg, cool_tconv, hicdc_anc_dict, w=5, graph=False, n=600, pc_rat = .1):
    s, e = map(int, [s, e])
    midpt = (chrom, s-5000*w, e+5000*w)
    right = (chrom, s, e+5000*(n-3))
    left = (chrom, s-5000*(n-3), e)

    expected_treg=treg_expecteds_balanced[chrom][0, 2:n]

    raw_L_reg = cool_treg.matrix(balance=True).fetch((midpt), (left))[:, :-2]
    raw_L_con = (cool_tconv.matrix(balance=True).fetch((midpt), (left)))[:, :-2]
    l = raw_L_reg.shape[1]
    mask = np.zeros(l)
    for i in range(l):
        anc = (chrom, s-i*5000, s-i*5000+5000)
        if hicdc_anc_dict.get(anc):
            mask[i-6:i+6] += 1
    mask = mask>4 
    mask_L = merge_mask(mask)[::-1]

    
    raw_R_reg = cool_treg.matrix(balance=True).fetch((midpt), (right))[:, 2:]
    raw_R_con = (cool_tconv.matrix(balance=True).fetch((midpt), (right)))[:, 2:]
    mask = np.zeros(l)
    for i in range(l):
        anc = (chrom, s+i*5000, s+i*5000+5000)
        if hicdc_anc_dict.get(anc):
            mask[i-6:i+6] += 1
    mask = mask>4    
    mask_R = merge_mask(mask)

    if graph==True:
        fig, axs = plt.subplots(1, 2)
        axs[0].plot(raw_L_reg.T)
        axs[0].set_title("Raw L reg")    

        axs[1].plot(raw_R_reg.T)
        axs[1].set_title("Raw R Reg")

        fig, axs = plt.subplots(1, 2)
        axs[0].plot(mask_L)
        axs[0].set_title("Mask L")    
        axs[1].plot(mask_R)
        axs[1].set_title("Mask R")
    
    
    row_L_reg = (do_interp(np.nanmean((raw_L_reg), axis=0),
                              )) + pc_rat*expected_treg[::-1]
    row_R_reg = (do_interp(np.nanmean((raw_R_reg), axis=0),
                              )) + pc_rat*expected_treg

    
    expected_tcon=tcon_expecteds_balanced[chrom][0, 2:n]
    row_L_con = (do_interp(np.nanmean(raw_L_con, axis=0),
                              )) + pc_rat*expected_tcon[::-1]
    
    row_R_con = (do_interp(np.nanmean(raw_R_con, axis=0),
                              )) + pc_rat*expected_tcon
    
    if graph==True:
        fig, axs = plt.subplots(1, 2)
        axs[0].plot(row_L_reg)
        axs[0].plot(row_L_con)
        axs[0].set_title("Row L")    
        axs[1].plot(row_R_reg)
        axs[1].plot(row_R_con)
        axs[1].set_title("Row R")
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
    
    
    expected_treg= expected_treg + pc_rat*expected_treg

    row_L_normed_reg = row_L_reg / expected_treg[::-1]
    row_R_normed_reg = row_R_reg / expected_treg
    row_L_normed_reg[np.isnan(row_L_normed_reg)] = 1
    row_R_normed_reg[np.isnan(row_R_normed_reg)] = 1

    
    expected_tcon= expected_tcon + pc_rat*expected_tcon

    row_L_normed_con = row_L_con / expected_tcon[::-1]
    row_R_normed_con = row_R_con / expected_tcon

    row_L_normed_con[np.isnan(row_L_normed_con)] = 1
    row_R_normed_con[np.isnan(row_R_normed_con)] = 1

    row_L_normed_con = sum_over_mask(mask_L, np.log2(row_L_normed_con))
    row_R_normed_con = sum_over_mask(mask_R, np.log2(row_R_normed_con))

    row_L_normed_reg = sum_over_mask(mask_L, np.log2(row_L_normed_reg))
    row_R_normed_reg = sum_over_mask(mask_R, np.log2(row_R_normed_reg))

    row_treg = np.append(row_L_normed_reg, row_R_normed_reg)    
    row_tcon = np.append(row_L_normed_con, row_R_normed_con)

    if graph==True:
        fig, axs = plt.subplots(1, 2)
        axs[0].plot(row_L_normed_reg)
        axs[0].plot(row_L_normed_con)
        axs[0].set_title("Row L")    
        axs[1].plot(row_R_normed_reg)
        axs[1].plot(row_R_normed_con)    
        axs[1].set_title("Row R")
    
    
    treg_v = row_treg
    tcon_v = row_tcon
    fc = treg_v-tcon_v
    fc[np.isinf(fc)] = 0
    fc[np.isnan(fc)] = 0 
    fc[np.isinf(fc)] = 0
    
    rs = {'expected_treg' : expected_treg,
    'expected_tcon' : expected_tcon,
    'raw_R_reg' : raw_R_reg,
    'raw_R_con' : raw_R_con,
    'row_L_reg' : row_L_reg,
    'row_L_con' : row_L_con,
    'row_R_reg' : row_R_reg,
    'row_R_con' : row_R_con,
    'row_R_normed_reg' : row_R_normed_reg,
    'row_R_normed_con' : row_R_normed_con,
    'row_L_normed_reg' : row_L_normed_reg,
    'row_L_normed_con' : row_L_normed_con,

    'row_treg' : row_treg,
    'row_tcon' : row_tcon,}
    
    return rs, fc


# def make_expected(significant_tads, cool_treg, cool_tconv, pseudocount=1, sigma=3, balance=False):
#     treg_list_of_expected = []
#     tconv_list_of_expected = []
#     print("Pseudocount: ", pseudocount)
#     print("Sigma: ", sigma)
    
#     max_size = 10000
#     import time
#     if balance==False:
#         for tad_ind in range(0, len(significant_tads)):
#         #     t1 = time.time()
#             chrom, s, e = significant_tads[int(tad_ind)][:3]
#             s, e = map(int, [s, e])
#             matrix = cool_treg.matrix(balance=False).fetch((chrom, s, e))
#             matrix[np.isnan(matrix)] = 0
#             matrix += pseudocount
#             matrix = matrix.astype(float)
#             mat1_smooth = scipy.ndimage.gaussian_filter(matrix, sigma=sigma)
#             size = mat1_smooth.shape[0]
#             expected = np.zeros(max_size)
#             for i in range(0, size):
#                 vals = np.diag(mat1_smooth, k=i)
#                 expected[i] = np.mean(vals)
#             treg_list_of_expected.append(isotonify(expected))


#             matrix = cool_tconv.matrix(balance=False).fetch((chrom, s, e))
#             matrix[np.isnan(matrix)] = 0

#             matrix += pseudocount
#             matrix = matrix.astype(float)
#             mat1_smooth = scipy.ndimage.gaussian_filter(matrix, sigma=sigma)
#             size = mat1_smooth.shape[0]
#             expected = np.zeros(max_size)
#             for i in range(0, size):
#                 vals = np.diag(mat1_smooth, k=i)
#                 expected[i] = np.mean(vals)
#             tconv_list_of_expected.append(isotonify(expected))    
#         treg_list_of_expected = np.asarray(treg_list_of_expected)
#         tconv_list_of_expected = np.asarray(tconv_list_of_expected)

#     elif balance==True:
#         for tad_ind in range(0, len(significant_tads)):
#         #     t1 = time.time()
#             chrom, s, e = significant_tads[int(tad_ind)][:3]
#             s, e = map(int, [s, e])
#             matrix = cool_treg.matrix(balance=balance).fetch((chrom, s, e))
#             matrix[np.isnan(matrix)] = 0
#             matrix += pseudocount
#             mat1_smooth = scipy.ndimage.gaussian_filter(matrix, sigma=sigma)
#             size = mat1_smooth.shape[0]
#             expected = np.zeros(max_size)
#             for i in range(0, size):
#                 vals = np.diag(mat1_smooth, k=i)
#                 expected[i] = np.mean(vals)
#             treg_list_of_expected.append(expected)

#             matrix = cool_tconv.matrix(balance=balance).fetch((chrom, s, e))
#             matrix[np.isnan(matrix)] = 0

#             matrix += pseudocount
#             mat1_smooth = scipy.ndimage.gaussian_filter(matrix, sigma=sigma)
#             size = mat1_smooth.shape[0]
#             expected = np.zeros(max_size)
#             for i in range(0, size):
#                 vals = np.diag(mat1_smooth, k=i)
#                 expected[i] = np.mean(vals)
#             tconv_list_of_expected.append(expected)    
#         treg_list_of_expected = np.asarray(treg_list_of_expected)
#         tconv_list_of_expected = np.asarray(tconv_list_of_expected)
        
#     return treg_list_of_expected, tconv_list_of_expected


def make_expected(df, tad1, tad2, chromsizes, res=5000, balanced=True):
    chrom = tad1[0]
    chromsize = chromsizes[chrom]
    n_bins = chromsize//res

    _, s1, e1 = tad1
    _, s2, e2 = tad2
    assert s1 < s2
    L_size = (e1-s1)//res
    R_size = (e2-s2)//res

    expected = np.zeros((L_size, R_size))
    test = np.zeros((expected.shape))

    averages = df[(df.region == chrom)]
 
    diag = (s2 - e1)//res
    bottom = diag
    diags = np.flip(np.indices(expected.shape)[0]) + (np.indices(expected.shape)[1]) + bottom
    n_bins_diags = n_bins - diags
    for val in np.unique(diags):
        if balanced==True:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['balanced.sum'])
        elif balanced==False:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['count.sum'])
    return expected/n_bins_diags

def make_expected_for_single_tad(df, tad1, chromsizes, res=5000, balanced=False):
    chrom = tad1[0]
    chromsize = chromsizes[chrom]
    n_bins = chromsize//res

    _, s1, e1 = tad1
    L_size = (e1-s1)//res

    expected = np.zeros((L_size, L_size))
    test = np.zeros((expected.shape))

    averages = df[(df.region == chrom)]

    diag = 0
    diags = np.triu(np.flip(np.indices(expected.shape)[0]) + (np.indices(expected.shape)[1]))
    diags += (diags).T - np.diag(np.diag(diags))-diags[0, 0]
    n_bins_diags = n_bins - diags
    for val in np.unique(diags):
        if balanced==True:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['balanced.sum'])
        elif balanced==False:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['count.sum'])
    return expected/n_bins_diags



def make_full_expected(df, tad1, tad2, res=5000, balanced=True):
    chrom = tad1[0]
    _, s1, e1 = tad1
    _, s2, e2 = tad2
    assert s1 < s2
    L_size = (e2-s1)//res
    R_size = (e2-s1)//res
    expected = np.zeros((L_size, R_size))
    test = np.zeros((expected.shape))

    averages = df[(df.region == chrom)]
 
    diag = (s2 - e1)//res
    bottom = diag
    diags = np.triu(np.flip(np.indices(expected.shape)[0]) + (np.indices(expected.shape)[1]))
    diags += (diags).T - np.diag(np.diag(diags))-diags[0, 0]

    for val in np.unique(diags):
        if balanced==True:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['balanced.avg'])
        elif balanced==False:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['count.sum'])
    return expected

def make_full_expected_avg(df, tad1, tad2, res=5000, balanced=True):
    chrom = tad1[0]
    _, s1, e1 = tad1
    _, s2, e2 = tad2
    assert s1 < s2
    L_size = (e2-s1)//res
    R_size = (e2-s1)//res
    expected = np.zeros((L_size, R_size))
    test = np.zeros((expected.shape))

    averages = df[(df.region == chrom)]
 
    diag = (s2 - e1)//res
    bottom = diag
    diags = np.triu(np.flip(np.indices(expected.shape)[0]) + (np.indices(expected.shape)[1]))
    diags += (diags).T - np.diag(np.diag(diags))-diags[0, 0]

    for val in np.unique(diags):
        if balanced==True:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['balanced.avg'])
        elif balanced==False:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['count.avg'])
    return expected






def isotonify(vector):
    n = len(vector)
    for i in range(n):
        val = vector[i]
        vector[i:][vector[i:] > val] = val
    return vector


import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

# import skimage
# from copy import deepcopy
# import pybedtools as pbt
# import itertools
# def merged_calling(tad, cool_treg, cool_tconv, treg_list_of_expected, tconv_list_of_expected, diff_tads, res=5000, graph=False, balance=False):
#     loops = {}
#     (chrom, start, e) = tad
#     start, e = int(start), int(e)
#     tadstring = f'{chrom} {start} {e}'
#     tad_as_pbt = pbt.BedTool(tadstring, from_string = True)
#     subtads = diff_tads.intersect(tad_as_pbt, u = True).sort()
#     n_subtads = len(subtads)

#     if balance==False:
#         for subtad in subtads:
#             chrom, sub_s, sub_e = subtad[:3]
#             sub_s, sub_e = int(sub_s), int(sub_e)
#             sub_s = sub_s - 5*res
#             sub_e = sub_e + 5*res
#             new_subtad = (chrom, sub_s, sub_e)
#             newloops, _, _ = differential_calling(new_subtad, cool_treg, cool_tconv, treg_list_of_expected, tconv_list_of_expected, res=5000, graph=graph, balance=balance)
#             pruned_loops = deepcopy(newloops)
#             for loop in newloops:
#                 l1, l2 = loop
#                 if l2[1] < l1[1]:
#                     del pruned_loops[loop]
#                     continue
#                 elif ((l1[1]-sub_s)//res < 2) or ((sub_e - l2[2])//res < 2):
#                     del pruned_loops[loop]
#                     continue
#             loops.update(pruned_loops)

#         if n_subtads > 1:

#             start = int(subtads[0][1])
#             end = int(subtads[n_subtads-1][2])
#             megatad = (chrom, start, end)
#             newloops, _, _ = differential_calling(megatad, cool_treg, cool_tconv, treg_list_of_expected, tconv_list_of_expected, res=5000, graph=graph, balance=balance)
            
#             pruned_loops = deepcopy(newloops)
#             for subtad in subtads:
#                 subtad_start = int(subtad[1])
#                 subtad_end = int(subtad[2])
#                 for loop in newloops:
#                     if not pruned_loops.get(loop):
#                         continue
#                     l1, l2 = loop
#                     if l2[1] < l1[1]:
#                         del pruned_loops[loop]
#                         continue                
#                     if ((l1[1] > subtad_end) and (l1[1] < subtad_end)):
#                         del pruned_loops[loop] 

#             loops.update(pruned_loops)
#     return loops

# from skimage.morphology import disk

# def differential_calling(tad, cool_treg, cool_tconv, treg_list_of_expected, tconv_list_of_expected, res=5000, graph=False, balance=False):
#     loops = {}
#     # print("Res is 5000")

#     (chrom, start, e) = tad
#     start, e = int(start), int(e)

#     if balance==False:
#         matrix_treg, smooth_treg,  = make_matrices(tad, cool_treg, pseudocount=3, sigma=1)
#         query_treg, dropoff_treg = make_query(smooth_treg)
#         query_treg = isotonify(query_treg)


#         neighbor_inds_treg = return_nearest_neighbors(treg_list_of_expected, query_treg)
#         expected_treg = calculate_expected(treg_list_of_expected, neighbor_inds_treg)
        
#         # print(np.min(expected_treg))
#         # print(np.min(query_treg), 'hi')

#         exp_treg_mat = key_to_matrix(expected_treg, len(matrix_treg))

#         matrix_tconv, smooth_tconv = make_matrices(tad, cool_tconv, pseudocount=3, sigma=1)
#         query_tconv, dropoff_tconv = make_query(smooth_tconv)
#         query_tconv = isotonify(query_tconv)

#         neighbor_inds_tconv = return_nearest_neighbors(tconv_list_of_expected, query_tconv)
#         expected_tconv = calculate_expected(tconv_list_of_expected, neighbor_inds_tconv)
#         exp_tconv_mat = key_to_matrix(expected_tconv, len(matrix_tconv))

        
#         treg_expr = smooth_treg / exp_treg_mat
#         tconv_expr = smooth_tconv / exp_tconv_mat


#         treg_mat = deepcopy(treg_expr).astype(float)
#         treg_mat = (treg_mat - np.min(treg_mat))
#         treg_mat = treg_mat/10
#         treg_mat[treg_mat > 1] = 1

#         tconv_mat = deepcopy(tconv_expr).astype(float)
#         tconv_mat = (tconv_mat - np.min(tconv_mat))
#         tconv_mat = tconv_mat/10
#         tconv_mat[tconv_mat > 1] = 1


#         filtered_treg = skimage.filters.rank.mean_bilateral(treg_mat, disk(10), s0=10, s1=10)
#         filtered_tconv = skimage.filters.rank.mean_bilateral(tconv_mat, disk(10), s0=10, s1=10)
        
#         m_treg = filtered_treg.astype(float) - filtered_tconv.astype(float)
#         # m_treg -= 1
#         # m_treg[m_treg < 0] = 0

#         blobs_treg = extrema.h_maxima(m_treg, 20)
#         blobs_treg = np.where(blobs_treg > 0)

#         m_tconv = -m_treg
#         # m_tconv -= 1
#         # m_tconv[m_tconv < 0] = 0

#         blobs_tconv = extrema.h_maxima(m_tconv, 20)
#         blobs_tconv = np.where(blobs_tconv > 0)


#         blobs_treg_alone = extrema.h_maxima(filtered_treg, 10)
#         blobs_treg_alone = np.where(blobs_treg_alone > 0)

#         blobs_tcon_alone = extrema.h_maxima(filtered_tconv, 10)
#         blobs_tcon_alone = np.where(blobs_tcon_alone > 0)

#         xs_treg = []
#         ys_treg = []

#         maxsize = matrix_tconv.shape[0]

#         for blob in zip(*blobs_treg):
#             y, x = blob
#             if len(blobs_treg_alone[0]) > 0:
#                 closest_peak_index = np.argmin(np.power((blobs_treg_alone[0]-y), 2) + np.power((blobs_treg_alone[1]-x), 2))
#                 closest_peak = blobs_treg_alone[0][closest_peak_index], blobs_treg_alone[1][closest_peak_index]
#                 if np.abs(y-closest_peak[0]) + np.abs(x-closest_peak[1]) < 20:
#                     y, x = closest_peak
#                     # print(np.abs(y-closest_peak[0]) + np.abs(x-closest_peak[1]))
#                     # print(blobs_treg_alone)
#                     pass
#                 else:
#                     continue

#             ind1, ind2 = int(x), int(y)

#             if (ind1 < 2) or (ind2 < 2) or (ind1 > maxsize-3) or (ind2 > maxsize-3):
#                 continue

#             if filtered_treg[ind1, ind2] < 10:
#                 continue                


#             if m_treg[ind1, ind2] < 15:
#                 continue                

#             # if (matrix_treg/matrix_tconv)[ind1-2:ind1+3, ind2-2:ind2+3].mean() < 1.2:
#                 # continue
#             xs_treg.append(x)
#             ys_treg.append(y)

#             genloc_s = start+y*res
#             genloc_e = start+x*res
#             l1 = (chrom, genloc_s, genloc_s+res)
#             l2 = (chrom, genloc_e, genloc_e+res)
#             loops[(l1, l2)] = np.log2((matrix_tconv)/(matrix_treg))[ind1, ind2]

#         xs_tcon = []
#         ys_tcon = []

#         for blob in zip(*blobs_tconv):
#             y, x = blob
#             if len(blobs_tcon_alone[0]) > 0:            
#                 closest_peak_index = np.argmin(np.power((blobs_tcon_alone[0]-y), 2) + np.power((blobs_tcon_alone[1]-x), 2))
#                 closest_peak = blobs_tcon_alone[0][closest_peak_index], blobs_tcon_alone[1][closest_peak_index]
#                 # print("poop:")
#                 # print(closest_peak)
#                 if np.abs(y-closest_peak[0]) + np.abs(x-closest_peak[1]) < 20:
#                     y, x = closest_peak
#                     # print(np.abs(y-closest_peak[0]) + np.abs(x-closest_peak[1]))
#                     # print(blobs_tcon_alone)
#                     pass
#                 else:
#                     continue
#             ind1, ind2 = int(x), int(y)

#             if (ind1 < 2) or (ind2 < 2) or (ind1 > maxsize-3) or (ind2 > maxsize-3):
#                 continue

#             if filtered_tconv[ind1, ind2] < 10:
#                 continue                


#             if m_tconv[ind1, ind2] < 15:
#                 continue
            
#             xs_tcon.append(x)
#             ys_tcon.append(y)
#             # print(x, y)
#             genloc_s = start+y*res
#             genloc_e = start+x*res
#             # print(genloc_s, genloc_e)
#             l1 = (chrom, genloc_s, genloc_s+res)
#             l2 = (chrom, genloc_e, genloc_e+res)
#             loops[(l1, l2)] = np.log2((matrix_tconv)/(matrix_treg))[ind1, ind2]


#     if graph==True:
#         fig, axs = plt.subplots(1, 3, figsize=(15, 5))
#         xs = np.arange(len(query_treg))
#         print(len(xs))
#         axs[0].scatter(xs, query_treg, label='Real')
#         # axs[0].scatter(xs, expected_treg[:len(query_treg)], label='Expected')

#         axs[0].legend()
#         axs[1].scatter(xs, expected_treg[:len(query_treg)], label='Difference')
#         axs[1].legend()

#         axs[2].scatter(xs, expected_treg[:len(query_treg)] / query_treg, label='Difference')
#         axs[2].legend()
#         axs[0].set_xlabel('Distance from diagonal')
#         axs[0].set_ylabel('Average contact frequency')
#         axs[1].set_xlabel('Distance from diagonal')
#         axs[1].set_ylabel('Average contact frequency')

#         axs[0].set_title('Contact frequency: real')
#         axs[1].set_title('Contact frequency: real - expected')

#         axs[0].set_ylim([0, 175])

#         plt.legend()

#         print("SHIIIIIT")

#         fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#         c = axs[0].matshow(treg_expr, cmap=cm.gist_heat_r, vmin=.5, vmax = 4)
#         axs[1].matshow(tconv_expr, cmap=cm.gist_heat_r, vmin=.5, vmax = 4)
#         fig.subplots_adjust(right=.8)
#         cbar_ax = fig.add_axes([.85, .15, .05, .7])
#         fig.colorbar(c, cax=cbar_ax)
#         fig.suptitle('Divide by expected')

#         fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#         c = axs[1].matshow(filtered_treg, cmap=cm.gist_heat_r, vmin=0, vmax = 260)
#         axs[0].matshow(filtered_tconv, cmap=cm.gist_heat_r, vmin=0, vmax = 260)
#         fig.subplots_adjust(right=.8)
#         cbar_ax = fig.add_axes([.85, .15, .05, .7])
#         fig.colorbar(c, cax=cbar_ax)
#         fig.suptitle('Subtract Conditions')

#         fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#         c = axs[1].matshow(m_tconv, cmap=cm.gist_heat_r, vmin=-50, vmax = 260)
#         axs[0].matshow(m_treg, cmap=cm.gist_heat_r, vmin=-50, vmax = 260)
#         fig.subplots_adjust(right=.8)
#         cbar_ax = fig.add_axes([.85, .15, .05, .7])
#         fig.colorbar(c, cax=cbar_ax)
#         fig.suptitle('Subtract Conditions')

#         fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#         axs[0].matshow(matrix_treg, cmap=cm.gist_heat_r, norm = colors.LogNorm(vmin=np.min(smooth_treg[smooth_treg > 0]), vmax=np.max(matrix_treg)))
#         axs[1].matshow(matrix_tconv, cmap=cm.gist_heat_r, norm = colors.LogNorm(vmin=np.min(smooth_treg[smooth_treg > 0]), vmax=np.max(matrix_treg)))
        
#         fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#         axs[0].matshow(matrix_treg, cmap=cm.gist_heat_r, norm = colors.LogNorm(vmin=np.min(smooth_treg[smooth_treg > 0]), vmax=np.max(matrix_treg)))
#         axs[1].matshow(matrix_tconv, cmap=cm.gist_heat_r, norm = colors.LogNorm(vmin=np.min(smooth_treg[smooth_treg > 0]), vmax=np.max(matrix_treg)))
        
#         for x, y in zip(xs_treg, ys_treg):
#             c = plt.Circle((x, y), 4, color='green', linewidth=2, fill=False)
#             axs[0].add_patch(c)

#         for x, y in zip(xs_tcon, ys_tcon):
#             c = plt.Circle((x, y), 4, color='green', linewidth=2, fill=False)
#             axs[1].add_patch(c)

#     # return loops, matrix_treg, matrix_tconv
#     return loops, filtered_treg, filtered_tconv



# def make_matrices(example, cool,  pseudocount=3, sigma=3, balance=False):
#     # print("Pseudocount: ", pseudocount)
    
#     # print("Sigma: ", sigma)
#     chrom, s, e = example
#     matrix = cool.matrix(balance=balance).fetch((chrom, s, e))
#     matrix[np.isnan(matrix)] = 0
#     matrix += pseudocount
#     matrix = matrix.astype(float)
#     mat1_smooth = scipy.ndimage.gaussian_filter(matrix, sigma=sigma)
#     return matrix, mat1_smooth

# def make_query(mat1_smooth):
#     assert mat1_smooth.shape[0] == mat1_smooth.shape[1] 
#     query_size = len(mat1_smooth)
#     # print("Query size: ", query_size)

#     query = np.zeros(query_size)
#     for i in range(0, query_size):
#         vals = np.diag(mat1_smooth, k=i)
#         query[i] = np.mean(vals)

#     full_size = mat1_smooth.shape[0]
#     full_dropoff = np.zeros(full_size)
#     for i in range(0, full_size):
#         vals = np.diag(mat1_smooth, k=i)
#         full_dropoff[i] = np.mean(vals)

#     return query, full_dropoff



# def make_matrix_query(example, cool,  pseudocount=3, sigma=3, balance=False):
#     print("Pseudocount: ", pseudocount)
    
#     print("Sigma: ", sigma)
#     chrom, s, e = example
#     matrix = cool.matrix(balance=balance).fetch((chrom, s, e))
#     matrix[np.isnan(matrix)] = 0

#     query_size = len(matrix)
#     print("Query size: ", query_size)

#     matrix += pseudocount
#     mat1_smooth = scipy.ndimage.gaussian_filter(matrix, sigma=sigma)

#     query = np.zeros(query_size)
#     for i in range(0, query_size):
#         vals = np.diag(mat1_smooth, k=i)
#         query[i] = np.mean(vals)

#     full_size = matrix.shape[0]
#     full_dropoff = np.zeros(full_size)
#     for i in range(0, full_size):
#         vals = np.diag(mat1_smooth, k=i)
#         full_dropoff[i] = np.mean(vals)

#     return matrix, mat1_smooth, query, full_dropoff

# def return_nearest_neighbors(list_of_expected, query):
#     query_size = len(query)
#     dists = np.sum(np.power(np.asarray(list_of_expected)[:, :query_size] - query, 2), axis=1)
#     pinds = np.argsort(dists)
#     return pinds


# def calculate_expected(list_of_expected, pinds, n_of_keys=40):
#     vals = list_of_expected[pinds[1:n_of_keys], :]
#     vals[vals == 0] = np.nan
#     expected = np.nanmean(vals, axis=0)
#     expected[np.isnan(expected)] = np.min(expected[~np.isnan(expected)])
#     return expected

# def key_to_matrix(key, size):
#     res = np.zeros([size]*2)
#     for i in range(size):
#         rows = np.arange(0, size-i)
#         cols = rows + i
#         res[rows, cols] = key[i]
#     return res + res.T - np.diag(np.diag(res))
    






