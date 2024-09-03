import time
import numpy as np
import scipy
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
import matplotlib.pyplot as plt
import bbi
import cooler
from copy import deepcopy
import scipy.cluster

def make_order2(matrix):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method='average', metric='cosine')
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)

        ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
        order = scipy.cluster.hierarchy.leaves_list(ordering)
    else:
        order = np.arange(len(matrix))
    # order = dendro['leaves']
    return order    

def make_nice(mat):
    mat = deepcopy(mat)
    mat[np.isnan(mat)]=0
    mat[np.isinf(mat)] = 0
    mat[:, 0] += 1e-10
    return mat

def make_expected_for_zscore(raw_mat):
    means = []
    for i in range(len(raw_mat)):
        m = np.nanmean(np.diag(raw_mat, k=i))
        if np.isnan(m):
            break
        means.append(m)
    meanmat = np.zeros(raw_mat.shape)
    diags = np.triu(np.flip(np.indices(raw_mat.shape)[0]) + (np.indices(raw_mat.shape)[1]))
    diags += (diags).T - np.diag(np.diag(diags))-diags[0, 0]
    exp = deepcopy(diags).astype(float)
    exp[np.isnan(raw_mat)] = np.nan
    for c, val in enumerate(means):
        exp[diags==c] = val
    stds = []
    for i in range(len(raw_mat)):
        m = np.nanstd(np.diag(raw_mat, k=i))
        if np.isnan(m):
            break
        if i == 0:
            stds.append(m)
        elif i > 300:
            stds.append(stds[-1])
        else:
            if m < stds[-1]:
                stds.append(m)
            else:
                stds.append(stds[-1])

    diags = np.triu(np.flip(np.indices(raw_mat.shape)[0]) + (np.indices(raw_mat.shape)[1]))
    diags += (diags).T - np.diag(np.diag(diags))-diags[0, 0]
    std_mat = deepcopy(diags).astype(float)
    std_mat[np.isnan(std_mat)] = np.nan
    for c, val in enumerate(stds):
        std_mat[diags==c] = val        
    return exp, std_mat

def make_diag_zscore(mat):
    mat[np.isinf(mat)] = 0 
#     mat[np.isnan(mat)] = 0
    exp, std = make_expected_for_zscore(mat)
    return (mat - 0)/std

def make_trans_zscore(mat):
    mat[np.isinf(mat)] = 0 
#     mat[np.isnan(mat)] = 0
    exp, std = np.nanmean(mat), np.nanstd(mat)
    return (mat -exp)/std

def make_full_zscore_mat(cool1, cool2, balance=False):
    n = sum([len(region_to_inds[x]) for x in parsed_chroms])
    mat = np.zeros((n, n))
    tot1 = 0
    for chrom1 in parsed_chroms:
        tot2 = 0

        for chrom2 in parsed_chroms:
            if chrom1 == chrom2:
                submat = cool1.matrix(balance=balance).fetch(chrom1, chrom2)
                submat2 = cool2.matrix(balance=balance).fetch(chrom1, chrom2)                
                diff_mat = submat-submat2
                f = make_diag_zscore(diff_mat)
                n1, n2 = f.shape
                mat[tot1:tot1+n1, tot2:tot2+n2] = f
            else:
                submat = cool1.matrix(balance=balance).fetch(chrom1, chrom2)
                submat2 = cool2.matrix(balance=balance).fetch(chrom1, chrom2)                
                diff_mat = submat-submat2
                f = make_trans_zscore(diff_mat)
                n1, n2 = f.shape
                mat[tot1:tot1+n1, tot2:tot2+n2] = f        
            tot2 += n2

        tot1 += n1
    return mat


path_to_cool = '../loop_plots/'

print("Importing coolers")
cool_treg_big = cooler.Cooler(path_to_cool + 'Treg_all.mcool::resolutions/250000')
cool_tconv_big = cooler.Cooler(path_to_cool + 'Tconv_all.mcool::resolutions/250000')

chromsizes = {}
for chrom, size in cool_tconv_big.chromsizes.iteritems():
    chromsizes[chrom] = size

print("Making parsed chroms")
parsed_chroms = [x for x in chromsizes if (x != 'M') and (x != 'Y')]


print("Making inds/regions")
inds_to_region = {}
region_to_inds = {}
for chrom in parsed_chroms:
    df = cool_treg_big.bins()[:]
    subdf = df[df['chrom'] == chrom]
    regions = list(zip(subdf['chrom'].values, map(int, subdf['start'].values), map(int, subdf['end'].values)))
    inds = np.arange(len(regions))
    inds_to_region[chrom] = regions
    region_to_inds[chrom] = dict(zip(regions, inds))
    
all_region_to_ind = {}
tot = 0
for chrom in parsed_chroms:
    df = cool_treg_big.bins()[:]
    subdf = df[df['chrom'] == chrom]
    regions = list(zip(subdf['chrom'].values, map(int, subdf['start'].values), map(int, subdf['end'].values)))
    inds = np.arange(len(regions)) + tot
    for r, i in zip(regions, inds):
        all_region_to_ind[r] = i
    tot += len(regions)


print("Making full zscore mat")
full_zscore_mat = make_full_zscore_mat(cool_treg_big, cool_tconv_big, balance=True)

print("Making order2")
order = make_order2(make_nice(full_zscore_mat))
import pickle
print("Pickle dumping")
pickle.dump(order, open( "order.p", "wb" ) )
print("Making picture")
fig, ax = plt.subplots(figsize=(12, 12))
ax.matshow(full_zscore_mat[order, :][:, order], cmap=cm.bwr, vmin=-4, vmax=4)
fig.savefig('./plots/ordered_zscore_mat.png')

print("Done with all")
