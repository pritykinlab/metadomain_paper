import pandas as pd
import numpy as np
import pickle

def str2inds(s):
    i1, i2 = s.split('_')
    return int(i1), int(i2)

arr = np.asarray
def inds_from_index(inds):
    return arr(list(pd.Series(inds).apply(str2inds)))

def loops_to_anchors(loops):
    anchors = []
    for i in loops:
        anchors.append(i[:3])
        anchors.append(i[3:6])
    return anchors
    
def make_int(l):
    return l[0], int(l[1]), int(l[2])

arr = np.asarray
def extend_s(s, d = 10_000):
    s = arr(s).astype(int)
    s -= d
    return s

def extend_e(e, d = 10_000):
    e = arr(e).astype(int)
    e += d
    return e

def extend_s_e(s, e, d = 10_000):
    s = extend_s(s, d = d)
    e = extend_e(e, d = d)
    return s, e


def make_loop_int(l):
    return l[0], int(l[1]), int(l[2]), l[3], int(l[4]), int(l[5])

def make_loop_str(l):
    return l[0], str(l[1]), str(l[2]), l[3], str(l[4]), str(l[5])

def make_str(l):
    return l[0], str(l[1]), str(l[2])

def slice_mat(mat, i1, i2, d):
    return mat[i1-d:i1+d+1, i2-d:i2+d+1]

def extend_l(l, d):
    return l[0], l[1]-d, l[2]+d

def add_chr(l):
    l = list(l)
    if 'chr' not in l[0]:
        l[0] = 'chr' + l[0]
        return tuple(l)
    else:
        return tuple(l)
    
def remove_chr(l):
    l = list(l)
    if 'chr' in l[0]:
        l[0] = l[0][3:]
        return tuple(l)
    else:
        return tuple(l)
        
def remove_chr_bedtool(bedtool):
    new_bedtool_list = []
    for i in bedtool:
        z = list(i) 
        chrom = z[0]
        if 'chr' == chrom[:3]:
            newchrom = chrom[3:]
        else:
            newchrom = chrom
        z[0] = newchrom
        new_bedtool_list.append(z)
    new_bedtool = pbt.BedTool(new_bedtool_list)
    return new_bedtool

def remove_chr_anc(anc):
    anc = list(anc)
    chrom = anc[0]
    if 'chr' == chrom[:3]:
        newchrom = chrom[3:]
    else:
        newchrom = chrom
    anc[0] = newchrom
    return anc

def remove_chr_list(bedtool):
    new_bedtool_list = []
    for i in bedtool:
        z = list(i) 
        chrom = z[0]
        if 'chr' == chrom[:3]:
            newchrom = chrom[3:]
        else:
            newchrom = chrom
        z[0] = newchrom
        new_bedtool_list.append(z)
    return new_bedtool_list

def notexpanded_tuple_to_grange(tup):
    chrom, s, e = tup
    return f'{chrom}:{s}-{e}'


def get_name(xind, ind_to_gene, default="None"):
    name = ind_to_gene.get(xind, [])
    name = [x for x in name if 'Gm' not in x]
    if len(name)==0:
        name=default
    else:
        name = name[0]
    return name

def format_pvalue(p, pco = .05, scico = .01,):
    if p == 0:
        p = '<1.8e-308'
    elif p > pco:
        p = 'NS'
    elif p > scico:
        p = str(np.round(p, 2))
    else:
        p = f'{p:.0e}'
    return p

def format_pval_as_asterisks(p, pco = .05, scico = .01, nsstr='NS'):
    if (p > .05) or (np.isnan(p)):
        p = nsstr
    elif p > .01:
        p = '*'
    elif p > .001:
        p = '**'
    else:
        p = '***'
    return p

def tuple_to_grange(chrom, s, e):
    return f'{chrom}:{s}-{e}'

def grange_to_tuple(grange):
    chrom, rest = grange.split(":")
    s, e = rest.split("-")
    return chrom, s, e

def unzreg(places):
    return list(zip(*places))

def get_unique(loops):
    return list(set(loops))

def add_chr_to_list(chroms):
    return ['chr'+ x for x in chroms]


def get_col(bed, col):
    return arr(list(zip(*bed))[col])

def get_bedtool_lengths(bedtool):
    es = get_col(bedtool, 2).astype(int) 
    ss = get_col(bedtool, 1).astype(int)
    deltas = np.abs(es-ss)
    return deltas

import scipy
from scipy.cluster.hierarchy import fcluster

def make_order_and_cluster_custom(matrix, method='average', metric='cosine', n_clusters=2):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric=metric)
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)
        order = dendro['leaves']
        print(n_clusters)
        clusters = scipy.cluster.hierarchy.fcluster(linkage, t=n_clusters, criterion='maxclust')-1
        return order, clusters, linkage
    else:
        order = np.arange(len(matrix))
        return [0]*len(order), order, None
    
    
def make_order_and_cluster_optimal(matrix, n_clusters, method='ward', metric='euclidean',
                                   return_linkage=False):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric=metric)
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)
        ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
        order = scipy.cluster.hierarchy.leaves_list(ordering)
        p = scipy.cluster.hierarchy.fcluster(linkage, t=n_clusters, criterion='maxclust')-1
    else:
        order = np.arange(len(matrix))
    # order = dendro['leaves']
    if return_linkage:
        return order, p, linkage
    else:
        return order, p


def overlap_two_separate_dfs(df1, df2):
    bed1, bed2 = index_to_bedtool(df1.index), index_to_bedtool(df2.index)
    l1s = []
    l2s = []
    for i in bed1.intersect(bed2, wo=True):
        l1 = i[:3]
        l2 = i[3:6]
        l1s.append(l1)
        l2s.append(l2)
    return bedtool_to_index(l1s), bedtool_to_index(l2s)

def rename_clusts_by_order(clusts, o):
    newclusts = clusts.copy()
    seen = set()
    for counter, i in enumerate(clusts[o]):
        if i in seen:
            continue
        else:
            seen.add(i)
            newclusts[clusts==i] = len(seen)-1
    return newclusts


def rename_clusts_by_clusts_in_order(clusts, clusts_in_order):
    newclusts = clusts.copy()
    for c, i in enumerate(clusts_in_order):
        newclusts[clusts==i] = c
    return newclusts


def filter_loops_by_dist(looplist, L=-np.inf, R=np.inf):
    filtloops = []
    for i in looplist:
        s1, s2 = i[1], i[4]
        s1, s2 = map(int, [s1, s2])
        if (abs(s2-s1) > L) and (abs(s2-s1) < R):
            filtloops.append(i)
        else:
            continue
    return set(filtloops)


def loops_intersecting_anchor(loops, anchor):
    l1s, l2s = [], []
    for i in loops:
        l1, l2 = i[:3], i[3:6]
        l1s.append(l1)
        l2s.append(l2)

    l1s_with_anchor = pbt.BedTool(l1s).intersect([anchor], c=True)
    l2s_with_anchor = pbt.BedTool(l2s).intersect([anchor], c=True)
    l1s_in_anchor = get_col(l1s_with_anchor, -1).astype(int)
    l2s_in_anchor = get_col(l2s_with_anchor, -1).astype(int)
    loops_with_anchor = []
    for i in range(len(loops)):
        if (l1s_in_anchor[i] > 0) | (l2s_in_anchor[i] > 0):
            loops_with_anchor.append(loops[i])
    return loops_with_anchor

import pybedtools as pbt
def remove_chr_bedtool_loops(bedtool):
    new_bedtool_list = []
    for i in bedtool:
        z = list(i) 
        chrom = z[0]
        if 'chr' == chrom[:3]:
            newchrom = chrom[3:]
        else:
            newchrom = chrom
        z[0] = newchrom

        chrom = z[3]
        if 'chr' == chrom[:3]:
            newchrom = chrom[3:]
        else:
            newchrom = chrom
        z[3] = newchrom
        new_bedtool_list.append(z)
    new_bedtool = pbt.BedTool(new_bedtool_list)
    return new_bedtool

def add_chr_to_anc(anc):
    anc = list(anc)
    if 'chr' not in anc[0]:
        anc[0] = 'chr' + anc[0]
    return tuple(anc)

def add_chr_to_bedtool(bedtool):
    newls = []
    for i in bedtool:
        l = list(i)
        if 'chr' not in l[0]:
            l[0] = 'chr' + l[0]
        newls.append(l)
    return pbt.BedTool(newls)

def make_ind_conversion(all_ind_to_region, all_ind_to_region_50kb):
    ind_to_depth = {}; ind_to_ind = {}; n = len(all_ind_to_region_50kb)
    for c, reg_250 in enumerate(all_ind_to_region):
        k = c*5
        U = min(n-1, k+200)
        L = max(0, k-200)
        for j in range(L, U):
            reg_50 = all_ind_to_region_50kb[j]
            is_contained = check_contained(reg_250, reg_50)
            if is_contained:
                start_250 = reg_250[1]
                start_50 = reg_50[1]
                depth = (start_50-start_250)//50_000
                ind_to_ind[j] = c
                ind_to_depth[j] = depth
                
    ind_to_ind_reverse = dict(zip(ind_to_ind.values(), ind_to_ind.keys()))
    return ind_to_depth, ind_to_ind, ind_to_ind_reverse

def check_contained(big_reg, small_reg):
    chrom1, s1, e1 = big_reg
    chrom2, s2, e2 = small_reg
    if (chrom1==chrom2) and (s2 >= s1) and (e2 <= e1):
        return True
    else:
        return False
    pass

def bedprop(a, b, frac=True):
    z = a.intersect(b, u=True)
    if frac:
        return len(z)/len(a)
    else:
        return len(z)

def transpose_dict(x):
    a, b = zip(*x.items())
    return dict(zip(b, a))
    
def bedcount(a, b):
    z = a.intersect(b, u=True)
    z2 = a.subtract(b, A=True)
    assert len(z) + len(z2) == len(a)
    return len(z), len(z2)

def loop_grange_to_granges(loop_grange):
    return loop_grange.split('|')

def loop_grange_to_tuple(loop_grange):
    grange1, grange2 = loop_grange.split('|')
    l1, l2 = map(grange_to_tuple, [grange1, grange2])
    return tuple(list(l1) + list(l2))


def granges_to_tuple_set(granges):
    return set(map(grange_to_tuple, granges))

def bedtool_to_grange_set(bedtool):
    grange_set = set()
    for i in bedtool:
        grange_set.add(tuple_to_grange(*i[:3]))
    return grange_set

def add_chr_to_loop_bedtool(loops):
    ls = []
    loops = list(loops)
    for i in loops:
        l1 = add_chr_to_anc(i[:3])
        l2 = add_chr_to_anc(i[3:6])
        l = list(l1) + list(l2)
        ls.append(l)
    return pbt.BedTool(ls)

def relabel_clusters(cluster, order):
    seen = set()
    newcluster = cluster.copy()
    for c, i in enumerate(cluster[order]):
        if i in seen:
            continue
        else:
            seen.add(i)
            newcluster[cluster==i] = len(seen)-1
    return newcluster

def loop_grange_list_to_anc_granges(loop_grange_list):
    anc_granges = set()
    for i in loop_grange_list:
        l1, l2 = loop_grange_to_granges(i)
        l1, l2 = map(grange_to_tuple, [l1, l2])
        anc_granges.add(l1)
        anc_granges.add(l2)
    return list(anc_granges)

def df_to_loop_grange_list(df):
    loop_grange_list = []
    for _, row in  df.iterrows():
        loop_grange = granges_to_loop_grange(row[:3], row[3:6])
        loop_grange_list.append(loop_grange)
    return loop_grange_list

def granges_to_loop_grange(grange1, grange2):
    return ('|').join([grange1, grange2])


def intersects_diagonal(l1, l2):
    if l2[1] < l1[2]:
        return True
    else:
        return False

def chromify_for_cool(cool, grange):
    if 'chr' in cool.chromnames[0]:
        return add_chr_to_anc(grange)
    else:
        return remove_chr(grange)

def get_differential_sites(df, pco=.05, lfc=0):
    sites_up = df.index[(df['padj'] < pco) & (df['log2FoldChange'] > lfc)]
    sites_down = df.index[(df['padj'] < pco) & (df['log2FoldChange'] < -lfc)]
    return sites_up, sites_down

def get_differential_sites_separate_df(lfc_df, pval_df, column, pco=.05, lfc=0):
    sites_up = lfc_df.index[(pval_df[column] < pco) & (lfc_df[column] > lfc)]
    sites_down = lfc_df.index[(pval_df[column] < pco) & (lfc_df[column] < -lfc)]
    return sites_up, sites_down

def get_pileup_from_bigwig(bigwig, places, delta=5000, bins=100, genome='mm10'):
    places = places.slop(b=delta, genome=genome)
    chrom, s, e = get_col(places, 0), get_col(places, 1).astype(int), get_col(places, 2).astype(int)
    return bigwig.stackup(chrom, s, e, bins=bins)

def get_mean_val_from_bigwig(bigwig, places, bins=1):
    chrom, s, e = get_col(places, 0), get_col(places, 1).astype(int), get_col(places, 2).astype(int)
    return bigwig.stackup(chrom, s, e, bins=bins)


def get_merged_mats(cooldict, chrom1, s1, e1, chrom2, s2, e2, delta):
    vdict = {}
    for cond in cooldict:
        m = cooldict[cond].matrix(balance=True).fetch((chrom1, s1-delta, e1+delta), (chrom2, s2-delta, e2+delta))
        vdict[cond] = m
    return vdict

def subset_dict(d, keys):
    newd = {}
    for k in keys:
        newd[k] = d[k]
    return newd

def in_chunks(data, n):
    """Yield successive n-sized chunks from data."""
    for i in range(0, len(data), n):
        yield data[i:i + n]

def load_pickle(pkl):
    with open(pkl, 'rb') as f:
        val = pickle.load(f)
    return val

def get_label_from_vs(vs):
    base, lfc = vs.split("_thresh=")
    base = ' ÷ '.join(base.split("_vs_")[::-1])
    return f'{base} ({lfc})\n'

def get_label_from_multiindex(col):
    return col[0].split("_")[0] + " " + col[1].split("_")[-1]

def fetch_mean_matched_values(needle, haystack, with_duplicates=False, n=1):
    inds = []
    for _, i in enumerate(needle):
        t = (haystack - i).abs()
        if n == 1:
            ord_inds = np.argmin(t)
            index_to_add = t.index[ord_inds]
            if with_duplicates:
                inds.append(index_to_add)
            else:
                if index_to_add not in inds:
                    inds.append(index_to_add)
        else:
            ord_inds = np.argsort(t)[:n]
            index_to_add = t.index[ord_inds]
            if with_duplicates:
                inds.extend(index_to_add)
            else:
                for i in index_to_add:
                    if i not in inds:
                        inds.append(i)
    matched_means = haystack.loc[list(inds)]
    return matched_means

def fetch_mean_matched_values_biased_up(basemean_series, series_to_match, with_duplicates=False):
    inds = []
    for _, i in enumerate(basemean_series):
        t = (series_to_match - i)
        t = t[t > 0].abs()
        if len(t) == 0:
            continue
        ord_inds = np.argmin(t)
        index_to_add = t.index[ord_inds]
        if with_duplicates:
            inds.append(index_to_add)
        else:
            if index_to_add not in inds:
                inds.append(index_to_add)
    matched_means = series_to_match.loc[list(inds)]
    return matched_means

def index_to_bedtool(index):
    return pbt.BedTool([x.split("_") for x in index])

def bedtool_to_index(bedtool, n=3):
    return pd.Series(['_'.join(x[:n]) for x in list(bedtool)])

def trans(x):
    return list(zip(*x))

def nonan_test(v1, v2, test=scipy.stats.ranksums):
    v1 = v1[~np.isnan(v1)]
    v2 = v2[~np.isnan(v2)]
    return test(v1, v2)

def samesize_nonan_test(v1, v2, test=scipy.stats.pearsonr):
    bad = np.isnan(v1) | np.isnan(v2) | np.isinf(v1) | np.isinf(v2)
    v1 = v1[~bad]
    v2 = v2[~bad]
    if (len(v1) < 4) or (len(v2) < 4):
        return np.nan, np.nan
    else:
        return test(v1, v2)    

from collections import defaultdict
def reverse_dict(original_dict):
    reversed_dict = defaultdict(list)
    for key, value in original_dict.items():
        if isinstance(value, list):
            for v in value:
                reversed_dict[v].append(key)
        else:
            reversed_dict[value].append(key)
    return reversed_dict

def df_to_bedtool(df):
    return pbt.BedTool.from_dataframe(df)


