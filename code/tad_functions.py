import numpy as np
import scipy
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
import matplotlib.pyplot as plt
import bbi

def make_mats(df, chrom, tad_index_mapping):
    indices = df['index'].values
    n = len(tad_index_mapping[chrom])
    lfcs = np.zeros((n, n))*np.nan
    pvals = np.zeros((n, n))*np.nan

    pvalue = df.pvalue.values
    lfc = df.log2FoldChange
    i1s, i2s = [], []
    for i in indices:
        i1, i2 = tad_index_mapping[chrom][i[0]], tad_index_mapping[chrom][i[1]]
        i1s.append(i1)
        i2s.append(i2)
    pvals[i1s, i2s] = pvalue
    pvals[i2s, i1s] = pvalue

    lfcs[i1s, i2s] = lfc
    lfcs[i2s, i1s] = lfc
    return pvals, lfcs


def make_int(tad):
    return (tad[0], int(tad[1]), int(tad[2]))

def make_str(tad):
    return (tad[0], str(tad[1]), str(tad[2]))


def make_cross_mats(df, chrom1, chrom2, tad_index_mapping):
    indices = df['index'].values
    n1 = len(tad_index_mapping[chrom1])
    n2 = len(tad_index_mapping[chrom2])
    lfcs = np.zeros((n1, n2))*np.nan
    pvals = np.zeros((n1, n2))*np.nan

    pvalue = df.pvalue.values
    lfc = df.log2FoldChange
    i1s, i2s = [], []
    for i in indices:
        i1, i2 = tad_index_mapping[chrom1][i[0]], tad_index_mapping[chrom2][i[1]]
        i1s.append(i1)
        i2s.append(i2)
    pvals[i1s, i2s] = pvalue
    lfcs[i1s, i2s] = lfc
    return pvals, lfcs

from copy import deepcopy
def make_clusters(lfcs):
    w = deepcopy(lfcs)
    w[np.isnan(w)] = 0
    z = (w - np.nanmean(w, axis=0))
    z /= np.sqrt(np.sum((z*z), axis=0))
    q = (z.T)@(z)

    # Make clustering
    clusters = np.zeros(q.shape)
    cluster = 0  
    cluster_list = []
    clustered = []
    for i in range(1, len(q)):
        if q[i-1, i] > .7:
            clustered.append(i-1)
        else:
            clustered.append(i-1)
            cluster_list.append(clustered)
            clustered = []
            if i == len(q)-1:
                clustered.append(i)
                cluster_list.append(clustered)
    for cluster, i in enumerate(cluster_list):
        inds = np.asarray(i)
        v1 = inds[0]
        v2 = inds[-1]
        clusters[v1:v2+1, v1:v2+1] = cluster
    return clusters

def score_matrix(lfcs, pvals, methods, method_name, kind='intra'):
    scores = []
    for i in range(pvals.shape[0]):
        pval_row = pvals[i, :]
        lfc_row = lfcs[i, :]
        lfc_row = lfc_row[~np.isnan(pval_row)]    
        pval_row = pval_row[~np.isnan(pval_row)]
        score = methods[method_name](pval_row, lfc_row)
        if kind == 'inter':
            if len(lfc_row) / len(lfcs[i, :]) < .3:
                score = 1
        scores.append(score)
    return np.asarray(scores)


def make_tracks(tad_to_index_merged, res_dict):
    tracks = {}
    xs = []
    ys = []
    for ind, i in enumerate(tad_to_index_merged):
        tad = i[0], str(i[1]), str(i[2])
        tad_lfcs = res_dict['tad_to_lfcs'][tad]
        for lfc in tad_lfcs:
            xs.append(ind)
            ys.append(lfc)
    tracks['gene_LFC'] = [xs, ys]

    xs = []
    ys = []
    for ind, i in enumerate(tad_to_index_merged):
        tad = i[0], str(i[1]), str(i[2])
        tad_lfc = np.mean(res_dict['tad_27ac'][tad])
        xs.append(ind)
        ys.append(tad_lfc)
    tracks['tad_27ac'] = [xs, ys]


    xs = []
    ys = []
    for ind, i in enumerate(tad_to_index_merged):
        tad = i[0], str(i[1]), str(i[2])
        tad_lfc = np.mean(res_dict['tad_atac'][tad])
        xs.append(ind)
        ys.append(tad_lfc)
    tracks['tad_LFC'] = [xs, ys]

    xs = []
    ys = []
    for ind, i in enumerate(tad_to_index_merged):
        tad = i[0], str(i[1]), str(i[2])
        tad_value = res_dict['tad_foxp3'][tad]
        xs.append(ind)
        ys.append(tad_value)
    tracks['tad_foxp3'] = [xs, ys]

    xs = []
    ys = []
    for ind, i in enumerate(tad_to_index_merged):
        tad = i[0], str(i[1]), str(i[2])
        tad_value = res_dict['tad_ctcf'][tad]
        xs.append(ind)
        ys.append(tad_value)
    tracks['tad_ctcf'] = [xs, ys]

    xs = []
    ys = []
    for ind, i in enumerate(tad_to_index_merged):
        tad = i[0], str(i[1]), str(i[2])
        tad_value = res_dict['tad_tcf1'][tad]
        xs.append(ind)
        ys.append(tad_value)
    tracks['tad_tcf1'] = [xs, ys]    
    return tracks


def make_foxp3_track(foxp3_peaks_yuri, atac, tad_to_index_merged):
    xs = []
    ys = []
    for ind, i in enumerate(tad_to_index_merged):
        tad = i[0], str(i[1]), str(i[2])
        tad_value = res_dict['tad_foxp3'][tad]
        xs.append(ind)
        ys.append(tad_value)
    tracks['tad_foxp3'] = [xs, ys]
    
    tad_foxp3 = {}
    for i in tadset.intersect(foxp3_peaks_yuri, c=True):
        tad = tuple(i[:3])
        dist = (int(tad[2])-int(tad[1]))//1000
        c = int(i[3])
        tad_foxp3[tad] = c/dist*100

def make_pearson(mat):
    pearson_chrom_1 = mat
    pearson_chrom_1[np.isnan(pearson_chrom_1)] = 0
    v = (pearson_chrom_1.T - np.mean(pearson_chrom_1, axis=1))
    v /= np.sqrt(np.sum(v*v, axis=0))
    pearson_1 = (v.T)@(v)
    return pearson_1





def make_track_input(tadset, peak_genes, foxp3_peaks_yuri, tcf1_rest, ctcf_peaks, tcon_v_treg, diff_27_ac):
    tad_to_lfcs = {}
    for i in tadset:
        tad = tuple(i[:3])
        tad_to_lfcs.setdefault(tad, [])
    for i in tadset.intersect(peak_genes, wo=True):
        tad = tuple(i[:3])
        pval = float(i[7])
        assert pval/np.abs(pval) != -1
        lfc = float(i[9]) 
        if pval < .05:
            tad_to_lfcs[tad].append(lfc)

        
    tad_tcf1 = {}
    for i in tadset.intersect(tcf1_rest, c=True):
        tad = tuple(i[:3])
        dist = (int(tad[2])-int(tad[1]))//1000
        c = int(i[3])
        tad_tcf1[tad] = c/dist*100    
        
        
    tad_ctcf = {}
    for i in tadset.intersect(ctcf_peaks, c=True):
        tad = tuple(i[:3])
        dist = (int(tad[2])-int(tad[1]))//1000
        c = int(i[3])
        tad_ctcf[tad] = c/dist*100    
 
    tad_foxp3 = {}
    for i in tadset.intersect(foxp3_peaks_yuri, c=True):
        tad = tuple(i[:3])
        dist = (int(tad[2])-int(tad[1]))//1000
        c = int(i[3])
        tad_foxp3[tad] = c/dist*100    


    tad_atac = {}
    for i in tadset.intersect(tcon_v_treg, wo=True):
        tad = tuple(i[:3])
        lfc = float(i[7])
        tad_atac.setdefault(tad, [])
        tad_atac[tad].append(lfc)

    tad_27ac = {}
    for i in tadset.intersect(diff_27_ac, wo=True):
        tad = tuple(i[:3])
        lfc = float(i[7])
        tad_27ac.setdefault(tad, [])
        tad_27ac[tad].append(lfc)

    for tad in tadset:
        tad = tuple(tad[:3])
        if not tad_27ac.get(tad):
            tad_27ac[tad] = 0 


    for tad in tadset:
        tad = tuple(tad[:3])
        if not tad_atac.get(tad):
            tad_atac[tad] = 0 

    res_dict = {
        'tad_to_lfcs' : tad_to_lfcs,
        'tad_atac' : tad_atac,
        'tad_foxp3' : tad_foxp3,
        'tad_ctcf' : tad_ctcf,
        'tad_tcf1' : tad_tcf1,
        'tad_27ac' : tad_27ac,
    }
    return res_dict


def reorder_tracks(tracks, order_dict):
    newtracks = {}
    for name in tracks:
        xs, ys = tracks[name]
        newxs = list(map(lambda x: [order_dict[x]], xs))
        newtracks[name] = [newxs, ys]
    return newtracks




def make_interaction_matrix(all_lfc_values, all_raw_values, significant_tads, all_tad_to_index,):
    n = len(significant_tads)
    interaction_matrix = np.zeros((n, n))
    for c1, i in enumerate(significant_tads):
        for c2, j in enumerate(significant_tads):
            chrom1, chrom2 = i[0], j[0]
            tad1 = i[0], int(i[1]), int(i[2])
            tad2 = j[0], int(j[1]), int(j[2])
            ind1, ind2 = all_tad_to_index[tad1], all_tad_to_index[tad2]
            lfc = all_lfc_values[ind1, ind2]

            interaction_matrix[c1, c2] = lfc
    interaction_matrix[np.isnan(interaction_matrix)] = 0 


    n = len(significant_tads)
    raw_matrix = np.zeros((n, n))
    for c1, i in enumerate(significant_tads):
        for c2, j in enumerate(significant_tads):
            chrom1, chrom2 = i[0], j[0]
            tad1 = make_int(i)
            tad2 = make_int(j)
            ind1, ind2 = all_tad_to_index[tad1], all_tad_to_index[tad2]
            lfc = all_raw_values[ind1, ind2]
            raw_matrix[c1, c2] = lfc
    raw_matrix[np.isnan(interaction_matrix)] = 0 
    return interaction_matrix, raw_matrix


def plot_interaction_matrix(interaction_matrix, raw_matrix, significant_tads, significant_interchrom_tads, 
                             chrom_colormap, binary_palette):
    chrom_colors = []
    for c1, i in enumerate(significant_tads):
        chrom = i[0]
        chrom_colors.append(chrom_colormap[chrom])
    chrom_colors = np.asarray(chrom_colors)

    inter_sig_colors = []
    for c1, i in enumerate(significant_tads):
        if i in significant_interchrom_tads:
            inter_sig_colors.append(binary_palette[1])
        else:
            inter_sig_colors.append(binary_palette[0])
    inter_sig_colors = np.asarray(inter_sig_colors)

    # both_matrix = np.concatenate([interaction_matrix, raw_matrix], axis=1)
    order = make_order(interaction_matrix)

    fig, ax = plt.subplots(1, 6, figsize=(24, 10), gridspec_kw={'width_ratios': [14, .5, .5, 14, .5, .5]})
    ax[0].matshow(interaction_matrix[order, :][:, order], cmap = cm.bwr, vmin=-4, vmax=4, aspect='auto')
    ax[0].set_title("LFC")
    ax[1].matshow(np.expand_dims(inter_sig_colors[order], axis=1), aspect='auto')
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_title("Significantly changing")
    ax[2].matshow(np.expand_dims(chrom_colors[order], axis=1), aspect='auto')
    ax[2].set_xticks([])
    ax[2].set_yticks([])
    ax[2].set_title("Chromosomes")

    ax[0+3].matshow(raw_matrix[order, :][:, order], cmap = cm.bwr, vmin=0, vmax=2, aspect='auto')
    ax[0+3].set_title("LFC")
    ax[1+3].matshow(np.expand_dims(inter_sig_colors[order], axis=1), aspect='auto')
    ax[1+3].set_xticks([])
    ax[1+3].set_yticks([])
    ax[1+3].set_title("Significantly changing")
    ax[2+3].matshow(np.expand_dims(chrom_colors[order], axis=1), aspect='auto')
    ax[2+3].set_xticks([])
    ax[2+3].set_yticks([])
    ax[2+3].set_title("Chromosomes")    
    return order


import pybedtools as pbt
def make_tad_to_gene(peak_genes, significant_tads):
    tad_to_gene = {}
    for tad in significant_tads:
        try:
            genelist = zip(*list(zip(*peak_genes.intersect(pbt.BedTool("\t".join([str(x) for x in tad])
                                                              , from_string=True), u=True)))[3:5])
        except Exception as e:
            print(e)
            continue
        genes = []
        pval_list = []
        for i in genelist:
            genes.append(i)
            pval_list.append(float(i[1]))
        if len(pval_list) > 0:
            ind = np.argmin(pval_list)
            gene = genes[ind][0]
            tad_to_gene[tad] = gene
        else:
            tad_to_gene[tad] = "NA"    
    return tad_to_gene


# def make_order(matrix):
#     linkage = scipy.cluster.hierarchy.linkage(matrix, metric='cosine')
#     dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
#                                         color_threshold=-np.inf)
#     order = dendro['leaves']
#     return order

def make_order(matrix):
    linkage = scipy.cluster.hierarchy.linkage(matrix, method='ward', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)
    order = dendro['leaves']
    return order    

def make_order2(matrix):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method='ward', metric='euclidean')
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)

        ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
        order = scipy.cluster.hierarchy.leaves_list(ordering)
    else:
        order = np.arange(len(matrix))
    # order = dendro['leaves']
    return order    



def make_df(v1, v2, name1, name2):
    df = pd.DataFrame()
    values = list(v1) + list(v2)
    labels = [name1]*len(v1) + [name2]*len(v2)
    df['values']=values
    df['labels']=labels
    return df    

def make_df_from_dict(dic, **kwargs):
    df = pd.DataFrame()
    values = []
    labels = []
    for key in dic:
        values = values + dic[key]
        labels = labels + [key]*len(dic[key])
    df['values']=values
    df['labels']=labels
    return df    



def manual_eval(st):
    test1 = st.split('\'')
    chrom1 = test1[1]
    chrom2 = test1[3]
    test2 = st.split(',')
    s1 = int(test2[1])
    e1 = int(test2[2][:-1])
    s2 = int(test2[4])
    e2 = int(test2[-1][:-2])
    tad1 = (chrom1, s1, e1)
    tad2 = (chrom2, s2, e2)
    return ((tad1, tad2))

def make_raw(parsed_chroms, index_to_tad_merged, tad_to_index_merged, expected_dict, chromsizes, res):
    tcon_val_dict = {}
    treg_val_dict = {}

    for chrom in parsed_chroms:
        df = pd.read_csv(f'./deseq_tad_interactions/raw_signal/{chrom}_{chrom}_MERGED_TADS.csv', sep='\t')
        df['index'] = df['Unnamed: 0'].apply(manual_eval)
        df['Treg'] = df.iloc[:, ["Treg" in x for x in df.columns]].sum(axis=1)
        df['Tcon'] = df.iloc[:, ["Tcon" in x for x in df.columns]].sum(axis=1)
        indices = df['index'].values
        n = len(tad_to_index_merged[chrom])
        treg_vals = np.zeros((n, n))*np.nan
        tcon_vals = np.zeros((n, n))*np.nan

        i1s, i2s = [], []
        for i in indices:
            i1, i2 = tad_to_index_merged[chrom].get(i[0], -1), tad_to_index_merged[chrom].get(i[1], -1)
            if i2 == -1 or i1 == -1:
                print(i)
                continue
            i1s.append(i1)
            i2s.append(i2)
        treg_vals[i1s, i2s] = df['Treg']
        treg_vals[i2s, i1s] = df['Treg']

        tcon_vals[i1s, i2s] = df['Tcon']
        tcon_vals[i2s, i1s] = df['Tcon']


        chromsize = chromsizes[chrom]
        n_bins = chromsize//res

        df = expected_dict['treg_expected']
        xs = (df[df['region'] == chrom].diag)*res
        diags =  n_bins - (df[df['region'] == chrom].diag)
        ys = (df[df['region'] == chrom]['count.sum'])/diags



        dist_to_exp = dict(zip(xs, ys))

        n = len(tcon_vals)
        expected_vals = np.zeros(tcon_vals.shape)
        for i in range(n):
            for j in range(n):
                tad1 = index_to_tad_merged[chrom][i]
                tad2 = index_to_tad_merged[chrom][j]
                s1, e1, = int(tad1[1]), int(tad1[2])
                s2, e2, = int(tad2[1]), int(tad2[2])        
                midpt1 = (s1+e1)/2
                midpt2 = (s2+e2)/2     
                dist = abs(midpt2-midpt1)//5000*5000
                exp = dist_to_exp[dist]*(e1-s1)/5000*(e2-s2)/5000
                expected_vals[i, j] = exp

        tcon_vals = (tcon_vals/expected_vals)
        tcon_vals[np.isinf(tcon_vals)] = 0
        tcon_val_dict[chrom] = tcon_vals
        
        treg_vals = (treg_vals/expected_vals)
        treg_vals[np.isinf(treg_vals)] = 0
        treg_val_dict[chrom] = treg_vals
    print("NOT USING LOG2")
    return tcon_val_dict, treg_val_dict


from copy import deepcopy
def make_b_a_vectors(significant_tads, tad_to_index_merged, intra, cross, size, size_inter, 
                    compartment_vector_dict, compartment_interchrom_vector_dict):

    B_lfcs = []
    cross_B_lfcs = []
    A_lfcs = []
    cross_A_lfcs = []
    for i in significant_tads:
        i = make_int(i)
        chrom = i[0]
        ind = tad_to_index_merged[chrom][i]
        lfc_vector = deepcopy(intra[chrom][ind, :])
        inter_lfc_vector = deepcopy(cross[chrom][ind, :])

        # lfc_vector[np.isnan(lfc_vector)] = 0
        # inter_lfc_vector[np.isnan(inter_lfc_vector)] = 0

        c_vec = compartment_vector_dict[chrom]
        A_vec = c_vec*(c_vec > 0)
        B_vec = c_vec*(c_vec < 0)

        size_vector = size[chrom]
        size_inter_vector = size_inter[chrom]
        
        c_vec = compartment_interchrom_vector_dict[chrom]
        A_inter_vec = c_vec*(c_vec > 0)
        B_inter_vec = c_vec*(c_vec < 0)

        good = ~np.isnan(lfc_vector)
        B_lfc = np.average(lfc_vector[good], weights=(B_vec*size_vector)[good])
        A_lfc = np.average(lfc_vector[good], weights=(A_vec*size_vector)[good])

        good = ~np.isnan(inter_lfc_vector)
        B_inter_lfc = np.average(inter_lfc_vector[good], weights=(B_inter_vec*size_inter_vector)[good])
        A_inter_lfc = np.average(inter_lfc_vector[good], weights=(A_inter_vec*size_inter_vector)[good])

        
        B_lfcs.append(B_lfc)
        cross_B_lfcs.append(B_inter_lfc)    
        A_lfcs.append(A_lfc)
        cross_A_lfcs.append(A_inter_lfc)    
        
    B_lfcs = np.asarray(B_lfcs)
    A_lfcs = np.asarray(A_lfcs)
    cross_B_lfcs = np.asarray(cross_B_lfcs)
    cross_A_lfcs = np.asarray(cross_A_lfcs)
    return B_lfcs, A_lfcs, cross_B_lfcs, cross_A_lfcs


def make_bedtool(bedtool, significant_tads, treg_ac_bw, tcon_ac_bw, ):
    sizes = []
    peak_to_tad = {}
    for i in bedtool.intersect(pbt.BedTool(significant_tads), wo=True):
        peak = tuple(i[:3])
        tad = tuple(i[-4:-1])
        peak = "chr" + peak[0], peak[1], peak[2]
        peak_to_tad[peak]=tad
        size = int(peak[2]) - int(peak[1])
        sizes.append(size)
    sizes = np.asarray(sizes)
        
    chroms, ss, es = list(zip(*peak_to_tad.keys()))

    treg_ac_scores = np.log2(treg_ac_bw.stackup(chroms, ss, es, bins=1000))
    tn_ac_scores = np.log2(tcon_ac_bw.stackup(chroms, ss, es, bins=1000))
    treg_ac_scores[np.isinf(treg_ac_scores)] = 0
    tn_ac_scores[np.isinf(tn_ac_scores)] = 0
    treg_ac_scores = np.nanmean(treg_ac_scores, axis=1)[:, None]*sizes
    tn_ac_scores = np.nanmean(tn_ac_scores, axis=1)[:, None]*sizes
    tad_scores = {}
    for i in significant_tads:
        tad_scores.setdefault(i, [])    
    for i, peak in enumerate(peak_to_tad):
        tad_scores[peak_to_tad[peak]].append(treg_ac_scores[i]-tn_ac_scores[i])
    scores = []
    for i in significant_tads:
        final_val = np.nanmean(tad_scores[i])
        if np.isnan(final_val):
            final_val = 0
        scores.append(final_val)
    return np.asarray(scores)


def make_mat(order, significant_tads, tad_to_index_merged, cross_lfc_dict):
    inter_lfcs = []
    for c, i in enumerate(order):
        tad = significant_tads[i]
        chrom = tad[0]
        ind = tad_to_index_merged[chrom][make_int(tad)]
        row = cross_lfc_dict[chrom][ind, :]
        inter_lfcs.append(row)

    max_inter_size = 0
    for i in inter_lfcs:
        max_inter_size = max(max_inter_size, len(i))

    inter_mat = np.zeros((len(inter_lfcs), max_inter_size))*np.nan
    for c, i in enumerate(inter_lfcs):
        inter_mat[c, :len(i)] = i
    return inter_mat


def make_compartment_vector(chrom, tad_to_index_merged, comp_bbi):
    chroms, ss ,es, = [], [],[]
    for i in tad_to_index_merged[chrom]:
        chrom, s, e = i[:3]
        chroms.append(chrom)
        ss.append(s)
        es.append(e)
    compartment_scores = np.nanmean(comp_bbi.stackup(chroms, ss, es, bins=1000), axis=1)
    return compartment_scores


def make_intra_sizes(order, significant_tads, tad_to_index_merged):
    all_sizes = []
    for c, i in enumerate(order):
        tad = significant_tads[i]
        chrom = tad[0]
        sizes = np.asarray([(int(x[2])-int(x[1]))/1000 for x in tad_to_index_merged[chrom]])
        all_sizes.append(sizes/np.sum(sizes))

    max_n = 0
    for i in all_sizes:
        max_n = max(max_n, len(i))
    weight_mat = np.zeros((len(all_sizes), max_n))*np.nan
    for c, i in enumerate(all_sizes):
        weight_mat[c, :len(i)] = i
    return weight_mat

def make_chrom_intra_sizes(chrom, tad_to_index_merged):
    all_sizes = []
    for i in tad_to_index_merged[chrom]:
        tad = i
        chrom = tad[0]
        sizes = (int(tad[2])-int(tad[1]))
        all_sizes.append(sizes)
    all_sizes = np.asarray(all_sizes)
    return all_sizes/all_sizes.sum()

def make_chrom_inter_sizes(chrom1, parsed_chroms, tad_to_index_merged):
    all_sizes = []
    for chrom2 in parsed_chroms:
        sizes = []
        if 'Y' in chrom2 or 'M' in chrom2 or chrom1==chrom2:
            continue
        if chrom1 == chrom2:
            continue
        else:
            for tad in tad_to_index_merged[chrom2]:
                start = int(tad[1])
                end = int(tad[2])
                dist = end-start
                sizes.append(dist)
        all_sizes += sizes
    all_sizes = np.asarray(all_sizes)
    return all_sizes/all_sizes.sum()

def make_inter_sizes(order, significant_tads, all_tad_to_index):
    all_sizes = []
    for c, i in enumerate(order):
        tad = significant_tads[i]
        chrom = tad[0]
        sizes = []
        for i in all_tad_to_index:
            _, s, e = i
            if _ == chrom:
                continue
            else:
                dist = int(e)-int(s)
                diag = dist/1000
                sizes.append(diag)
        sizes = np.asarray(sizes)
        all_sizes.append(sizes/np.sum(sizes))
    max_n = 0
    for i in all_sizes:
        max_n = max(max_n, len(i))
    weight_mat = np.zeros((len(all_sizes), max_n))*np.nan
    for c, i in enumerate(all_sizes):
        weight_mat[c, :len(i)] = i
    return weight_mat


def make_self_raw(tad, expected, cooler):
    vals = cooler.matrix(balance=False).fetch(tad).astype(float)
    vals[np.isnan(expected)] = np.nan
    enrichment = np.nansum(vals)/np.nansum(expected)
    return enrichment

def make_self_vals(order, significant_tads, tad_to_index_merged, intra_lfcs):
    inter_lfcs = []
    for c, i in enumerate(order):
        tad = significant_tads[i]
        chrom = tad[0]
        ind = tad_to_index_merged[chrom][make_int(tad)]
        val = intra_lfcs[chrom][ind, ind]
        inter_lfcs.append(val)
    return np.expand_dims(inter_lfcs, axis=1)



def plot_region(cool, chrom, start, end, balance, norm):
    v = cool.matrix(balance=balance).fetch((chrom, start, end))
    fig, ax = plt.subplots(figsize=(8, 8))
    if norm == True:
        ax.matshow(v, cmap=cm.gist_heat_r, norm=colors.LogNorm())    
    else:
        ax.matshow(v, cmap=cm.gist_heat_r)
    return v


def get_fcs(regions, peak_genes):
    all_vals = list(zip(*peak_genes.intersect(pbt.BedTool(regions), u=True)))
    if len(all_vals) > 0:
        names = all_vals[3]    
        fcs = np.asarray(all_vals[6]).astype(float)
        pvals = np.asarray(all_vals[4]).astype(float)
        return names, fcs, pvals
    else: 
        return np.asarray(0), np.asarray(0), np.asarray(0)

def make_nice(mat):
    mat = deepcopy(mat)
    mat[np.isnan(mat)]=0
    mat[np.isinf(mat)] = 0
    mat[:, 0] += 1e-10
    return mat

def do_regional_roc_analysis(all_raw_values, query_tad_subset, retrieve_tad_subset, all_tad_to_index):

    merged_data = np.squeeze(np.nanmean(all_raw_values[query_tad_subset, :], axis=0), axis=0)
    merged_data[np.isnan(merged_data)] = 0
    print(merged_data.shape)        
    pearson_rs = []
    for i in all_raw_values:
        i[np.isnan(i)] = 0
        i[np.isinf(i)] = 0
        pearson_r = scipy.stats.pearsonr(i, merged_data)[0]
        pearson_rs.append(pearson_r)
    pearson_rs = np.asarray(pearson_rs)    

    scores = [0]
    tpr = []
    fpr = []
    sort_indices = np.argsort(-pearson_rs)
    n_tads = len(sort_indices)

    n = len(retrieve_tad_subset)
    for i in sort_indices:
        if i in retrieve_tad_subset:
            scores.append(scores[-1]+1)
        else:
            scores.append(scores[-1]+0)
        tpr.append(scores[-1]/n)
        fpr.append((len(scores)-scores[-1])/(n_tads-n))
    return scores, tpr, fpr, pearson_rs[sort_indices]



def do_roc_analysis(all_raw_values, significant_tads, query_tad_subset, retrieve_tad_order, all_tad_to_index):
    merged_data = np.zeros(len(all_raw_values))

    retrieve_tad_subset = [all_tad_to_index[make_int(significant_tads[x])] for x in retrieve_tad_order]
    for i in query_tad_subset:
        tad = significant_tads[i]
        ind = all_tad_to_index[make_int(tad)]
        vals = all_raw_values[ind, :]
        vals[np.isnan(vals)] = 0
        vals[np.isinf(vals)] = 0
        merged_data += vals
        
    pearson_rs = []
    for i in all_raw_values:
        i[np.isnan(i)] = 0
        i[np.isinf(i)] = 0
        pearson_r = scipy.stats.pearsonr(i, merged_data)[0]
        pearson_rs.append(pearson_r)
    pearson_rs = np.asarray(pearson_rs)    

    scores = [0]
    tpr = []
    fpr = []
    sort_indices = np.argsort(-pearson_rs)
    n_tads = len(sort_indices)

    n = len(retrieve_tad_subset)
    for i in sort_indices:
        if i in retrieve_tad_subset:
            scores.append(scores[-1]+1)
        else:
            scores.append(scores[-1]+0)
        tpr.append(scores[-1]/n)
        fpr.append((len(scores)-scores[-1])/(n_tads-n))
    return scores, tpr, fpr, pearson_rs[sort_indices]













def stouffer_naive(pval_row, lfc_row):
    pval_score = scipy.stats.combine_pvalues(pval_row, method='stouffer')[1]
    if pval_score == 0 :
        pval_score = 1e-300
    return np.log10(pval_score)

def count_significant(pval_row, lfc_row):
    count = np.mean(pval_row < .05)
    return count


def count_significant_lfc_no_pval(pval_row, lfc_row):
    count = np.mean((np.abs(lfc_row) > 2))
    return count

def count_significant_lfc(pval_row, lfc_row):
    count = np.mean((np.abs(lfc_row) > 2)& (pval_row < .05))
    return count

def fisher(pval_row, lfc_row):
    pval_score = scipy.stats.combine_pvalues(pval_row, method='fisher')[1]
    if pval_score == 0 :
        pval_score = 1e-300
    return np.log10(pval_score)

def stouffer_weight(pval_row, lfc_row):
    weights = np.abs(lfc_row)
    pval_score = scipy.stats.combine_pvalues(pval_row, method='fisher', weights=weights)[1]
    if pval_score == 0 :
        pval_score = 1e-300
    return np.log10(pval_score)

def stouffer_naive_remove_outliers(pval_row, lfc_row):
    pval_row[pval_row == 0] = 1e-300
    pval_row[pval_row == 1] = .9999
    pval_score = scipy.stats.combine_pvalues(pval_row, method='stouffer')[1]
    if pval_score == 0 :
        pval_score = 1e-300
    return np.log10(pval_score)

def fisher_remove_outliers(pval_row, lfc_row):
    pval_row[pval_row == 0] = 1e-300
    pval_row[pval_row == 1] = .9999
    pval_score = scipy.stats.combine_pvalues(pval_row, method='fisher')[1]
    if pval_score == 0 :
        pval_score = 1e-300
    return np.log10(pval_score)

def lfc_sum(pval_row, lfc_row):
    return np.nanmean(np.abs(lfc_row))

def stouffer_weight_remove_outliers(pval_row, lfc_row):
    pval_row[pval_row == 0] = 1e-300
    pval_row[pval_row == 1] = .9999
    weights = np.abs(lfc_row)
    pval_score = scipy.stats.combine_pvalues(pval_row, method='fisher', weights=weights)[1]
    if pval_score == 0 :
        pval_score = 1e-300
    return np.log10(pval_score)

def pval_test(lfc_row, **kwargs):
    pval_row = kwargs['pval_row']
    pval_row = pval_row[~np.isnan(pval_row)]
    return np.mean(pval_row < .05)

def lfc_pval_test_geom(lfc_row, **kwargs):
    lfc_row =lfc_row[~np.isnan(lfc_row)]
    sigma, mu = kwargs['sigma'], kwargs['mu']
    zscores = -np.abs((lfc_row-mu)/sigma)
    pvals_test = 2*scipy.stats.norm.cdf(zscores)
    mean_pval = np.exp(np.mean(np.log(pvals_test)))
    if len(zscores) < 40:
        mean_pval = 1
    return mean_pval

def lfc_pval_test(lfc_row, **kwargs):
    lfc_row =lfc_row[~np.isnan(lfc_row)]
    sigma, mu = kwargs['sigma'], kwargs['mu']
    zscores = -np.abs((lfc_row-mu)/sigma)
    pvals_test = 2*scipy.stats.norm.cdf(zscores)
    # mean_pval = np.exp(np.mean(np.log(pvals_test)))
    mean_pval = np.mean(pvals_test)
    if len(zscores) < 40:
        mean_pval = 1
    return mean_pval


def zscore_test(lfc_row, **kwargs):
    lfc_row =lfc_row[~np.isnan(lfc_row)]
    sigma, mu = kwargs['sigma'], kwargs['mu']
    zscores = np.abs((lfc_row-mu)/sigma)
    # pvals_test = 2*scipy.stats.norm.cdf(zscores)
    # mean_pval = np.exp(np.mean(np.log(zscores)))
    mean_pval = np.mean(zscores)
    if len(zscores) < 40:
        mean_pval = 0
    return mean_pval    

def empirical_test(lfc_row, **kwargs):
    lfc_row =lfc_row[~np.isnan(lfc_row)]
    interp_func = kwargs['interp_func']
    percentiles = interp_func(lfc_row)
    scores = 1-2*np.abs(percentiles-.5)
    mean_pval = np.exp(np.mean(np.log(scores)))
    if len(scores) < 40:
        mean_pval = 1
    return mean_pval    

# def empirical_pval_mean(lfc_row, **kwargs):
#     lfc_row =lfc_row[~np.isnan(lfc_row)]
#     interp_func = kwargs['interp_func']
#     percentiles = interp_func(lfc_row)
#     scores = 1-2*np.abs(percentiles-.5)
#     mean_pval = (scores < .1).mean()
#     return mean_pval    




def zscore_test_geom(lfc_row, **kwargs):
    lfc_row =lfc_row[~np.isnan(lfc_row)]
    sigma, mu = kwargs['sigma'], kwargs['mu']
    zscores = np.abs((lfc_row-mu)/sigma)
    # pvals_test = 2*scipy.stats.norm.cdf(zscores)
    mean_pval = np.exp(np.mean(np.log(zscores)))
    # mean_pval = np.mean(zscores)
    if len(zscores) < 40:
        mean_pval = 0
    return mean_pval    


def calculate_scores(lfcs, pvals, func, **kwargs,):
    scores = []
    n = len(lfcs)
    for i in range(lfcs.shape[0]):
        pval_row = pvals[i, :]
        lfc_row = lfcs[i, :]
        lfc_row = lfc_row[~np.isnan(lfc_row)]
        kwargs['pval_row'] = pval_row
        score = func(lfc_row, **kwargs)
        scores.append(score)
    scores = np.asarray(scores)
    return scores

def make_L(i, distances, min_tads=100, min_reach=10_000):
    temp = 0
    ind_counter = i
    while (temp < min_reach) and (ind_counter > 0):
        ind_counter -= 1
        temp += distances[ind_counter]

    if i - ind_counter < min_tads:
        ind_counter = max(0, i-min_tads)

    if ind_counter != 0:
        # print(i)
        # print(i-ind_counter)
        # print(np.sum(distances[ind_counter:i]))
        assert np.sum(distances[ind_counter:i]) >= min_reach and i-ind_counter  >= min_tads
    return ind_counter

def make_order2_and_cluster(matrix, n_clusters, method='ward', metric='euclidean'):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric=metric)
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)
        ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
        order = scipy.cluster.hierarchy.leaves_list(ordering)
        poop = fcluster(linkage, t=n_clusters, criterion='maxclust')-1
    else:
        order = np.arange(len(matrix))
    # order = dendro['leaves']
    return order, poop

def make_R(i, distances, min_tads=100, min_reach=10_000):
    n = len(distances)
    temp = 0
    ind_counter = i+1
    while (temp < min_reach) and (ind_counter < n):
        temp += distances[ind_counter]
        ind_counter += 1
    if ind_counter - i  < min_tads:
        ind_counter = min(n, i+min_tads)
    if ind_counter != n:
        # print(ind_counter-i)
        # print(np.sum(distances[i:ind_counter]))
        assert np.sum(distances[i:ind_counter]) > min_reach and ind_counter-i  >= min_tads
    return ind_counter    

    
def calculate_adjacent_scores(lfcs, pvals, func, distance=100, **kwargs,):
    scores = []
    n = len(lfcs)

    distances = kwargs['distances']
    for i in range(pvals.shape[0]):
        R = make_R(i, distances)
        L = make_L(i, distances)
        pval_row = pvals[i, L:R]
        lfc_row = lfcs[i, L:R]
        lfc_row = lfc_row[~np.isnan(lfc_row)]

        kwargs['pval_row'] = pval_row
        score = func(lfc_row, **kwargs)
        scores.append(score)
    scores = np.asarray(scores)
    return scores

def get_distances(tadlist):
    dists=[]
    for i in tadlist:
        d = (int(i[2])-int(i[1]))/5000
        dists.append(d)
    return np.asarray(dists)


def process_tads(merged_df, tad_to_index_merged, index_to_tad_merged, chromsizes, 
                func, func_name, threshold_func, threshold_dict, **kwargs):
    significant_tads = []
    intra_lfcs = {}
    intra_pvals = {}
    score_dict = {}
    threshold = threshold_dict[func_name]
    for chrom in chromsizes.keys():
        if 'Y' in chrom or 'M' in chrom:
            continue    
        distances = get_distances(index_to_tad_merged[chrom])
        kwargs['distances']=distances
        chrom_df = merged_df[merged_df['index'].apply(lambda x: x[0][0] == chrom)]
        pvals, lfcs = make_mats(chrom_df, chrom, tad_to_index_merged)

        intra_lfcs[chrom] = lfcs
        intra_pvals[chrom] = pvals

        # scores = calculate_adjacent_scores(lfcs, pvals, func, **kwargs)
        scores = calculate_scores(lfcs, pvals, func, **kwargs)
        significant = threshold_func(scores)

        if chrom =='4' or chrom=='1' or chrom == '8':
            fig, axs = plt.subplots(1, 2, figsize=(16, 12), gridspec_kw={'width_ratios': [3, 1]})
            axs[0].matshow(lfcs, cmap=cm.bwr, vmin=-2, vmax=2, aspect='auto',)
            xs = np.arange(len(scores))
            axs[1].plot(scores, xs)
            for a in axs[1:]:
                a.set_ylim([0, len(scores)])
                a.invert_yaxis()
                a.yaxis.set_major_locator(MaxNLocator(integer=True))
                a.vlines(threshold, 0, len(scores))
            fig.suptitle(chrom)
            fig.savefig(f'./plots/tad_interactions/chrom_plots/{chrom}_{func_name}.png')
            plt.close(fig)
        score_dict[chrom] = scores

        for j in np.where(significant)[0]:
            significant_tads.append(index_to_tad_merged[chrom][j])    

    return significant_tads, intra_lfcs, intra_pvals, score_dict


from scipy.interpolate import interp1d
def make_interp_func(lfcs):
    empiricals = []
    lfcs_for_empirical = np.ravel(lfcs[~np.isnan(lfcs)])
    L = np.min(lfcs_for_empirical)-1
    R = np.max(lfcs_for_empirical)+1
    for i in np.linspace(L, R, 1000):
        empirical = (lfcs_for_empirical < i).mean()
        empiricals.append(empirical)
    interp_func = interp1d(np.linspace(L, R, 1000), empiricals)
    return interp_func


def process_cross_tads(intra_lfcs, lfc_interchrom_dict, tad_to_index_merged, index_to_tad_merged, parsed_chroms,
        chromsizes, func, func_name, threshold_func, **kwargs):
    significant_interchrom_tads = []
    cross_lfc_dict = {}
    for chrom1 in parsed_chroms:        
        n1 = len(tad_to_index_merged[chrom1])
        cross_lfcs = np.zeros((n1, 0))*np.nan
        cross_pvals = np.zeros((n1, 0))*np.nan
        for chrom2 in parsed_chroms:
            if chrom1 < chrom2:
                lfc_mat = lfc_interchrom_dict[(chrom1, chrom2)]
            elif chrom1 == chrom2:
                continue
            else:
                lfc_mat = lfc_interchrom_dict[(chrom2, chrom1)].T
            cross_lfcs = np.concatenate((cross_lfcs, lfc_mat), axis=1)
        cross_lfc_dict[chrom1] = cross_lfcs
        scores = []
        for i in range(cross_pvals.shape[0]):
            lfc_row = cross_lfcs[i, :]
            lfc_row = lfc_row[~np.isnan(lfc_row)]    
            score = func(lfc_row, **kwargs)
            if len(lfc_row) / len(cross_lfcs[i, :]) < .1:
                score = func(lfc_row*0, **kwargs)
            scores.append(score)
        scores = np.asarray(scores)    

        significant = threshold_func(scores)
        for j in np.where(significant)[0]:
            significant_interchrom_tads.append(index_to_tad_merged[chrom1][j])

        if chrom1 == '1' or chrom1 == '4':
            fig, axs = plt.subplots(1, 3, figsize=(15, 6), gridspec_kw={'width_ratios': [4, 4, 2]})
            z = intra_lfcs[chrom1]
            axs[0].matshow(z, cmap=cm.bwr, vmin=-2, vmax=2, aspect='auto',)

            y = cross_lfcs
            axs[1].matshow(y, cmap=cm.bwr, vmin=-2, vmax=2, aspect='auto',)

            xs = np.arange(len(scores))
            axs[2].plot(scores, xs)
            for a in axs[2:]:
                a.set_ylim([0, len(scores)])
                a.invert_yaxis()
                a.yaxis.set_major_locator(MaxNLocator(integer=True))
            
            fig.suptitle(chrom1)
            fig.savefig(f'./plots/tad_interactions/chrom_plots/both_{chrom1}_{func_name}.png')
            plt.close(fig)
    return significant_interchrom_tads, cross_lfc_dict, 




import time
# def empirical_lfc(lfc_row, lfcs):
#     percentiles = []
#     for i in lfc_row:
#         if i < 0:
#             percentile = 2*((lfcs < i).mean())
#         if i > 0:
#             percentile = 2*((lfcs > i).mean())
#         percentiles.append(percentile)
        
#     mean_percentile = np.exp(np.mean(np.log(percentiles)))
#     # mean_pval = np.mean(pvals_test)
#     if len(percentiles) < 40:
#         mean_percentile = 1
#     return mean_percentile


methods = {    
    'lfc_sum' : lfc_sum,
    'count_significant' : count_significant,
    'count_significant_lfc' : count_significant_lfc,
    'fisher' : fisher,
    'fisher_remove_outliers' : fisher_remove_outliers,
    'stouffer_naive' : stouffer_naive,
    'stouffer_naive_remove_outliers' : stouffer_naive_remove_outliers,
    'stouffer_weight' : stouffer_weight,
    'stouffer_weight_remove_outliers' : stouffer_weight_remove_outliers,
    'count_significant_lfc_no_pval' : count_significant_lfc_no_pval,
    'lfc_pval_test' : lfc_pval_test,
}


def make_all_lfc_values(lfc_interchrom_dict, intra_lfcs, all_tad_to_index, tad_to_index_merged, ):
    n = len(all_tad_to_index)
    all_lfc_values = np.zeros((n, n))
    for tad1 in all_tad_to_index:
        for tad2 in all_tad_to_index:
            chrom1, chrom2 = tad1[0], tad2[0]
            i1, i2 = all_tad_to_index[tad1], all_tad_to_index[tad2]
            og_index1, og_index2 = tad_to_index_merged[chrom1][tad1], tad_to_index_merged[chrom2][tad2]

            if chrom1 < chrom2:
                lfc = lfc_interchrom_dict[(chrom1, chrom2)][og_index1, og_index2]
            elif chrom1 > chrom2:
                continue
            else:
                lfc = intra_lfcs[(chrom1)][og_index1, og_index2]
            all_lfc_values[i1, i2] = lfc
            all_lfc_values[i2, i1] = lfc
    return all_lfc_values



def make_all_raw_values(parsed_chroms, trans_exp_dict, tcon_val_dict, treg_val_dict, all_tad_to_index, tad_to_index_merged, 
                    chromsizes):
    n = len(all_tad_to_index)
    all_raw_values = np.zeros((n, n))
    all_raw_values_treg = np.zeros((n, n))

    raw_interchrom_dict = {}
    raw_interchrom_dict_treg = {}
    for chrom1 in parsed_chroms:
        for chrom2 in parsed_chroms:
            if chrom1 != chrom2 and chrom1 < chrom2:
                n1, n2 = map(len, [tad_to_index_merged[chrom1], tad_to_index_merged[chrom2]])
                raw_interchrom_dict[(chrom1, chrom2)] = np.zeros((n1, n2))
                raw_interchrom_dict_treg[(chrom1, chrom2)] = np.zeros((n1, n2))
                df = pd.read_csv(f'./deseq_tad_interactions/interchrom_raw/{chrom1}_{chrom2}.csv', sep='\t')
                tads = df['Unnamed: 0'].values
                tcon_vals = df.iloc[:, ["Tcon" in x for x in df.columns]].sum(axis=1).values+.5
                for i in range(len(df)):
                    tad1, tad2 = manual_eval(tads[i])
                    chrom1, chrom2 = tad1[0], tad2[0]
                    tot = trans_exp_dict[(chrom1, chrom2)]
                    s1, e1 = map(int, tad1[1:])
                    s2, e2 = map(int, tad2[1:])
                    size1, size2 = (e1-s1)//5000,  (e2-s2)//5000
                    total_size = size1*size2
                    both_chromsize = chromsizes[chrom1]/5000*chromsizes[chrom2]/5000
                    z = tot/both_chromsize*total_size
                    obs_div_exp = np.log2(tcon_vals[i]/(tot/both_chromsize*total_size))
                    i1, i2 = all_tad_to_index[tad1], all_tad_to_index[tad2]
                    all_raw_values[i1, i2] = obs_div_exp
                    all_raw_values[i2, i1] = obs_div_exp
                    i1_OG, i2_OG = tad_to_index_merged[chrom1][tad1], tad_to_index_merged[chrom2][tad2]
                    raw_interchrom_dict[(chrom1, chrom2)][i1_OG, i2_OG] = obs_div_exp
                raw_interchrom_dict[(chrom2, chrom1)] = raw_interchrom_dict[(chrom1, chrom2)].T

                treg_vals = df.iloc[:, ["Treg" in x for x in df.columns]].sum(axis=1).values+.5
                for i in range(len(df)):
                    tad1, tad2 = manual_eval(tads[i])
                    chrom1, chrom2 = tad1[0], tad2[0]
                    tot = trans_exp_dict[(chrom1, chrom2)]
                    s1, e1 = map(int, tad1[1:])
                    s2, e2 = map(int, tad2[1:])
                    size1, size2 = (e1-s1)//5000,  (e2-s2)//5000
                    total_size = size1*size2
                    both_chromsize = chromsizes[chrom1]/5000*chromsizes[chrom2]/5000
                    obs_div_exp = np.log2(treg_vals[i]/(tot/both_chromsize*total_size))
                    i1, i2 = all_tad_to_index[tad1], all_tad_to_index[tad2]
                    all_raw_values_treg[i1, i2] = obs_div_exp
                    all_raw_values_treg[i2, i1] = obs_div_exp

            if chrom1 == chrom2:
                tadlist = tad_to_index_merged[chrom1].keys()
                og_indices = [tad_to_index_merged[chrom1][make_int(tad)] for tad in tadlist]
                new_indices = [all_tad_to_index[make_int(tad)] for tad in tadlist]
                min_ind = np.min(new_indices)
                max_ind = np.max(new_indices)
                all_raw_values[min_ind:max_ind+1, min_ind:max_ind+1] = np.log2(tcon_val_dict[chrom1])
                all_raw_values_treg[min_ind:max_ind+1, min_ind:max_ind+1] = np.log2(treg_val_dict[chrom1])
    return raw_interchrom_dict, all_raw_values, all_raw_values_treg

def make_full_raw_interchrom_dict(parsed_chroms, tad_to_index_merged, raw_interchrom_dict):
    full_raw_interchrom_dict = {}
    for chrom1 in parsed_chroms:
        n1 = len(tad_to_index_merged[chrom1])
        temp_mat = np.zeros((n1, 0))
        for chrom2 in parsed_chroms:
            if chrom1 == chrom2:
                continue
            temp_mat = np.concatenate([temp_mat, raw_interchrom_dict[(chrom1, chrom2)]], axis=1)
        full_raw_interchrom_dict[chrom1] = temp_mat    
    return full_raw_interchrom_dict


def make_cluster_vals(clusters, matrix):
    us = np.unique(clusters)
    n = len(us)
    vals = np.zeros((n, n))
    for i in us:
        for j in us:
            rows = clusters==i
            cols = clusters==j
            vals[i, j] = np.nanmean(matrix[rows, :][:, cols])
    return vals

def make_reorder(clusters, matrix):
    us = np.unique(clusters)
    n = len(us)
    vals = np.zeros((n, n))
    
    for i in us:
        for j in us:
            rows = np.where(clusters==i)[0]
            cols = np.where(clusters==j)[0]
            vals[i, j] = (matrix[rows, :][:, cols]).sum()

    pearson = make_pearson(vals)
    new_order = make_order2(pearson)
    mapping = dict(zip(new_order, us))


    new_clusters = deepcopy(clusters)    
    for i in us:
        new_clusters[clusters==i] = mapping[i]
    return new_clusters

from scipy.cluster.hierarchy import fcluster
def make_order_and_cluster(mat, n_clusters):
    linkage = scipy.cluster.hierarchy.linkage(mat, method='ward', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)
    # ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, mat)
    # order = scipy.cluster.hierarchy.leaves_list(ordering)

    order = dendro['leaves']
    poop = fcluster(linkage, t=n_clusters, criterion='maxclust')-1

    poop = make_reorder(poop, make_pearson(mat))  
    return order, poop


from sklearn.cluster import KMeans
def make_order_and_kmeans_cluster(mat, n_clusters):
    linkage = scipy.cluster.hierarchy.linkage(mat, method='average', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)
    # ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, mat)
    # order = scipy.cluster.hierarchy.leaves_list(ordering)

    order = dendro['leaves']
    poop = KMeans(n_clusters=n_clusters).fit(mat).labels_

    poop = make_reorder(poop, mat)  
    return order, poop




def make_order_from_clusters(clusters, matrix):
    us = np.unique(clusters)
    n = len(clusters)
    tot = 0
    mapping = {}
    new_order = []
    for u in us:
        inds = clusters==u
        order = make_order2(matrix[inds, :])
        sub_ind_to_OG_ind = dict(zip(np.arange(inds.sum()), np.arange(len(inds))[inds]))
        new_order += [sub_ind_to_OG_ind[x] for x in order]
        tot += inds.sum()            
    return new_order


from matplotlib.ticker import MaxNLocator

def make_compartment_correlation(comp_vec, mat):
    rs = []
    for row in mat:
        r = scipy.stats.pearsonr(row, comp_vec)[0]
        rs.append((r))
    return np.asarray(rs)

def plot_raw_lfc_sidebyside(mat, raw_mat, order, chrom, compartment_vector_dict, treg_compartment_vector_dict, clusters=None):

    mat[np.isnan(mat)] = 0
    raw_mat[np.isinf(raw_mat)] = 0

    raw_mat = raw_mat[order, :][:, order]
    mat = mat[order, :][:, order]
    tcon_comp = compartment_vector_dict[chrom][order][None, :]
    treg_comp = treg_compartment_vector_dict[chrom][order][None, :]
    delta = (treg_comp - tcon_comp)

    raw_comp_corr = make_compartment_correlation(tcon_comp[0, :], raw_mat)
    lfc_comp_corr = make_compartment_correlation(tcon_comp[0, :], mat)


    fig, ax = plt.subplots(7, 2, figsize=(16, 8), gridspec_kw={'height_ratios': [20, 1, 1, 1, 1, 1, 1]} )    
    col = 0
    ax[0, col].matshow(mat, cmap=cm.bwr, vmin=-2, vmax=2, aspect='auto')
    ax[0, col].set_title(chrom)
    ax[1, col].matshow(tcon_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    ax[2, col].matshow(treg_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    if clusters is not None:

        k = np.max(clusters)
        clusters = clusters[order][None, :]
        ax[2, col].matshow(clusters, cmap=cm.tab20, vmin=0, vmax=k, aspect='auto')
    ax[3, col].matshow(delta, cmap=cm.bwr, vmin=-.01, vmax=.01, aspect='auto')
    ax[4, col].matshow(np.diag(mat)[None, :], cmap=cm.bwr, vmin=-.5, vmax=.5, aspect='auto')
    ax[5, col].matshow(raw_comp_corr[None, :], cmap=cm.bwr, vmin=-1, vmax=1, aspect='auto')
    ax[6, col].matshow(lfc_comp_corr[None, :], cmap=cm.bwr, vmin=-1, vmax=1, aspect='auto')    

    ylabels = ['Tn comp', 'Treg comp.', '∆', 'Self-LFC', 'Raw_comp_corr', 'LFC_comp_corr']
    for c, a in enumerate(ax[1:, col]):
        a.set_xticks([])
        a.set_yticks([])  
        a.set_ylabel(ylabels[c], rotation=0, loc='top')
        a.yaxis.set_label_position("right")
        a.yaxis.tick_right()

    col = 1
    ax[0, col].matshow(raw_mat, cmap=cm.bwr, vmin=-3, vmax=3, aspect='auto')
    ax[0, col].set_title(chrom)
    ax[1, col].matshow(tcon_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    ax[2, col].matshow(treg_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    if clusters is not None:
        ax[2, col].matshow(clusters, cmap=cm.tab20, vmin=0, vmax=k, aspect='auto')
    ax[3, col].matshow(delta, cmap=cm.bwr, vmin=-.01, vmax=.01, aspect='auto')
    ax[4, col].matshow(np.diag(mat)[None, :], cmap=cm.bwr, vmin=-.5, vmax=.5, aspect='auto')
    ax[5, col].matshow(raw_comp_corr[None, :], cmap=cm.bwr, vmin=-1, vmax=1, aspect='auto')
    ax[6, col].matshow(lfc_comp_corr[None, :], cmap=cm.bwr, vmin=-1, vmax=1, aspect='auto')    
    for a in ax[1:, col]:
        a.set_xticks([])
        a.set_yticks([])    
    return fig


def plot_all_clustered_raw(mat, raw_mat, order, chrom, all_compartment_vals, all_compartment_vals_treg, clusters=None):

    mat[np.isnan(mat)] = 0
    raw_mat[np.isinf(raw_mat)] = 0

    raw_mat = raw_mat[order, :][:, order]
    mat = mat[order, :][:, order]
    tcon_comp = all_compartment_vals[order][None, :]
    treg_comp = all_compartment_vals_treg[order][None, :]
    delta = (treg_comp - tcon_comp)

    k = np.max(clusters)
    clusters = clusters[order][None, :]

    fig, ax = plt.subplots(7, 2, figsize=(16, 8), gridspec_kw={'height_ratios': [20, 1, 1, 1, 1, 1, 1]} )    
    col = 0
    ax[0, col].matshow(mat, cmap=cm.bwr, vmin=-2, vmax=2, aspect='auto')
    ax[0, col].set_title(chrom)
    ax[1, col].matshow(tcon_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    ax[2, col].matshow(treg_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    ax[3, col].matshow(delta, cmap=cm.bwr, vmin=-.01, vmax=.01, aspect='auto')
    ax[4, col].matshow(clusters, cmap=cm.tab20, vmin=0, vmax=k, aspect='auto')

    ylabels = ['Tn comp', 'Treg comp.', '∆', 'Clusters', 'Raw_comp_corr', 'LFC_comp_corr']
    for c, a in enumerate(ax[1:, col]):
        a.set_xticks([])
        a.set_yticks([])  
        a.set_ylabel(ylabels[c], rotation=0, loc='top')
        a.yaxis.set_label_position("right")
        a.yaxis.tick_right()

    col = 1
    ax[0, col].matshow(raw_mat, cmap=cm.bwr, vmin=-3, vmax=3, aspect='auto')
    ax[0, col].set_title(chrom)
    ax[1, col].matshow(tcon_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    ax[2, col].matshow(treg_comp, cmap=cm.bwr, vmin=-.02, vmax=.02, aspect='auto')
    
    ax[3, col].matshow(delta, cmap=cm.bwr, vmin=-.01, vmax=.01, aspect='auto')
    ax[4, col].matshow(clusters, cmap=cm.tab20, vmin=0, vmax=k, aspect='auto')
    for a in ax[1:, col]:
        a.set_xticks([])
        a.set_yticks([])    
    return fig    

def plot_significant_by_significant(order, clusters, tad_to_gene, all_tad_to_index, significant_tads, tad_to_index_merged, tcon_val_dict,
        full_raw_interchrom_dict, self_raw_values_sig, tn_compartment_scores, treg_compartment_scores, intra_lfcs, 
        cross_lfc_dict, B_raw, A_raw, cross_B_raw, cross_A_raw, B_lfcs, A_lfcs, cross_B_lfcs, cross_A_lfcs, raw_matrix, interaction_matrix
    ):

    chroms = []
    for tad in significant_tads:
        chrom = tad[0]
        chroms.append(chrom)
    chroms = np.asarray(chroms)
    chroms = chroms[order]
    us  = np.unique(chroms)
    mapping = dict(zip(us, range(len(us))))
    new_chroms = [mapping[x] for x in chroms]
    chroms = np.asarray(new_chroms)

    intra_vals = make_mat(order, significant_tads, tad_to_index_merged, tcon_val_dict)
    weight_mat = make_intra_sizes(order, significant_tads, tad_to_index_merged)

    weight_mat_inter = make_inter_sizes(order, significant_tads, all_tad_to_index)
    inter_vals = make_mat(order, significant_tads, tad_to_index_merged, full_raw_interchrom_dict)
    inter_vals[np.isinf(inter_vals)] = 0

    self_vals = self_raw_values_sig[order][:, None]

    intra_raw_total = np.nansum(intra_vals*weight_mat, axis=1)[:, None]
    b_raw_vec = B_raw[order][:, None]
    a_raw_vec = A_raw[order][:, None]
    cross_b_raw_vec = cross_B_raw[order][:, None]
    cross_a_raw_vec = cross_A_raw[order][:, None]

    tn_vec = tn_compartment_scores[order]
    treg_vec = treg_compartment_scores[order]
    delta = (treg_compartment_scores-tn_compartment_scores)[order]
    sign_delta = (np.sign(treg_compartment_scores)-np.sign(tn_compartment_scores))[order]

    fig, ax = plt.subplots(1, 13, figsize=(17, 7), gridspec_kw={'width_ratios': [7, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5]} )

    sns.heatmap(raw_matrix[order, :][:, order], cmap=cm.bwr, vmin=-6, vmax=6, ax=ax[0], cbar=False,
                )
    sns.heatmap(chroms[:, None], cmap=cm.Accent, ax=ax[1], cbar=False)

    sns.heatmap(clusters[order][:, None], cmap=cm.tab20, ax=ax[2], cbar=False, vmin=0, vmax=10)
    sns.heatmap(intra_raw_total, cmap=cm.bwr, vmin=0, vmax=2, ax=ax[2+1], cbar=False)
    sns.heatmap(b_raw_vec, cmap=cm.bwr, vmin=0, vmax=2, ax=ax[3+1], cbar=False)
    sns.heatmap(cross_b_raw_vec, cmap=cm.bwr, vmin=0, vmax=2, ax=ax[4+1], cbar=False)
    sns.heatmap(a_raw_vec, cmap=cm.bwr, vmin=0, vmax=2, ax=ax[5+1], cbar=False)
    sns.heatmap(cross_a_raw_vec, cmap=cm.bwr, vmin=0, vmax=2, ax=ax[5+1+1], cbar=False)
    sns.heatmap(self_vals, vmin=0, vmax=2, cmap=cm.bwr, ax=ax[5+2+1], cbar=False)
    sns.heatmap(tn_vec, vmin=-4, vmax=4, cmap=cm.bwr, ax=ax[6+2+1], cbar=False)
    sns.heatmap(treg_vec, vmin=-4, vmax=4, cmap=cm.bwr, ax=ax[7+2+1], cbar=False)
    sns.heatmap(delta, vmin=-4, vmax=4, cmap=cm.bwr, ax=ax[8+2+1], cbar=False)
    sns.heatmap(sign_delta, vmin=-2, vmax=2, cmap=cm.bwr, ax=ax[8+2+1+1], cbar=False)



    ax[1+2].set_title("Intra")
    ax[1+3].set_title("Intra \n B raw")
    ax[1+4].set_title("Inter \n B raw")
    ax[1+5].set_title("Intra \n A raw")
    ax[1+6].set_title("Inter \n A raw")
    ax[1+5+2].set_title("Self")
    ax[1+6+2].set_title("Tn \n comp")
    ax[1+7+2].set_title("Treg \n comp")
    ax[1+8+2].set_title("Delta")
    ax[1+8+3].set_title("Sign ∂")    
    for a in ax[2:]:
        a.set_xticks([])
        a.set_yticks([])
    names = []
    for i in order:
        name = tad_to_gene[significant_tads[i]]
        names.append(name)

    ax[0].yaxis.set_major_locator(MaxNLocator(nbins=40, integer=True))    
    _ = ax[0].set_yticklabels(names)


    intra_vals = make_mat(order, significant_tads, tad_to_index_merged, intra_lfcs)    
    inter_vals = make_mat(order, significant_tads, tad_to_index_merged, cross_lfc_dict)
    self_vals = make_self_vals(order, significant_tads, tad_to_index_merged, intra_lfcs)

    intra_lfcs_total = np.nansum(intra_vals*weight_mat, axis=1)[:, None]
    b_lfc_vec = B_lfcs[order][:, None]
    a_lfc_vec = A_lfcs[order][:, None]
    cross_b_lfc_vec = cross_B_lfcs[order][:, None]
    cross_a_lfc_vec = cross_A_lfcs[order][:, None]

                
    fig, ax = plt.subplots(1, 13, figsize=(17, 7), gridspec_kw={'width_ratios': [7, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5]} )
    sns.heatmap(interaction_matrix[order, :][:, order], cmap=cm.bwr, vmin=-6, vmax=6, ax=ax[0], cbar=False)
    sns.heatmap(chroms[:, None], cmap=cm.Accent, ax=ax[1], cbar=False)
    sns.heatmap(clusters[order][:, None], cmap=cm.tab20, ax=ax[2], cbar=False, vmin=0, vmax=10)
    sns.heatmap(intra_lfcs_total, cmap=cm.bwr, vmin=-1, vmax=1, ax=ax[2+1], cbar=False)
    sns.heatmap(b_lfc_vec, cmap=cm.bwr, vmin=-2, vmax=2, ax=ax[3+1], cbar=False)
    sns.heatmap(cross_b_lfc_vec, cmap=cm.bwr, vmin=-2, vmax=2, ax=ax[4+1], cbar=False)
    sns.heatmap(a_lfc_vec, cmap=cm.bwr, vmin=-2, vmax=2, ax=ax[5+1], cbar=False)
    sns.heatmap(cross_a_lfc_vec, cmap=cm.bwr, vmin=-2, vmax=2, ax=ax[5+1+1], cbar=False)
    sns.heatmap(self_vals, vmin=-2, vmax=2, cmap=cm.bwr, ax=ax[5+ 2+1], cbar=False)
    sns.heatmap(tn_vec, vmin=-4, vmax=4, cmap=cm.bwr, ax=ax[6+2+1], cbar=False)
    sns.heatmap(treg_vec, vmin=-4, vmax=4, cmap=cm.bwr, ax=ax[7+2+1], cbar=False)
    sns.heatmap(delta, vmin=-4, vmax=4, cmap=cm.bwr, ax=ax[8+2+1], cbar=False)
    sns.heatmap(sign_delta, vmin=-2, vmax=2, cmap=cm.bwr, ax=ax[8+2+2], cbar=False)

    ax[1+2].set_title("Intra")
    ax[1+3].set_title("Intra \n B LFC")
    ax[1+4].set_title("Inter \n B LFC")
    ax[1+5].set_title("Intra \n A LFC")
    ax[1+6].set_title("Inter \n A LFC")
    ax[1+5+2].set_title("Self")
    ax[1+6+2].set_title("Tn \n comp")
    ax[1+7+2].set_title("Treg \n comp")
    ax[1+8+2].set_title("Delta")
    ax[1+8+3].set_title("∂(sign)")


    for a in ax[:]:
        a.set_xticks([])
        a.set_yticks([])
        
        
    names = []
    for i in order:
        name = tad_to_gene[significant_tads[i]]
        names.append(name)

    ax[0].yaxis.set_major_locator(MaxNLocator(nbins=40, integer=True))    
    _ = ax[0].set_yticklabels(names)



def chipseq_barplot(chrom, index_to_tad_merged, chrom_cluster_dict, bedtool_dict, atac_count,  all_tcon_v_treg):
    # Normalized by ATAC-seq
    val_dict = {}
    for name in bedtool_dict:
        bedtool = bedtool_dict[name]
        clusts = chrom_cluster_dict[chrom]
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
            
            for i in pbt.BedTool(changing_tads).intersect(bedtool.intersect(all_tcon_v_treg, u=True), c=True):
                tad = tuple(i[:3])
                count = int(i[3])
                if count > 0:
                    lfc_list.append(count/(atac_count[tad]))
                else:
                    lfc_list.append(0)
            val_dict[u] = lfc_list
        test_df = make_df_from_dict(val_dict)
        fig, ax = plt.subplots()
        sns.barplot(x="labels", y="values", data=test_df, ax=ax)
        ax.set_title(name)
        top = np.max([np.mean(x) for x in val_dict.values()])
        ax.set_ylim([0, 2*top])
        fig.savefig(f'./plots/tad_interactions/chrom_chipseq/{chrom}_{name}_atac.png')
        plt.close(fig)


    val_dict = {}
    for name in bedtool_dict:
        bedtool = bedtool_dict[name]
        clusts = chrom_cluster_dict[chrom]
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
            
            for i in pbt.BedTool(changing_tads).intersect(bedtool.intersect(all_tcon_v_treg, u=True), c=True):
                dist = (int(i[2])-int(i[1]))/1000                
                tad = tuple(i[:3])
                count = int(i[3])
                count = count/dist
                lfc_list.append(count)
            val_dict[u] = lfc_list
        test_df = make_df_from_dict(val_dict)

        fig, ax = plt.subplots()
        sns.barplot(x="labels", y="values", data=test_df, ax=ax)
        ax.set_title(name)
        top = np.max([np.mean(x) for x in val_dict.values()])
        ax.set_ylim([0, 2*top])
        fig.savefig(f'./plots/tad_interactions/chrom_chipseq/{chrom}_{name}.png')
        plt.close(fig)



def genes_barplot(chrom, chrom_cluster_dict, peak_genes, index_to_tad_merged):
    bedtool = peak_genes
    val_dict = {}
    clusts = chrom_cluster_dict[chrom]
    for u in np.unique(clusts):
        lfc_list = []
        inds = np.where(clusts == u)[0]
        changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
        for i in pbt.BedTool(changing_tads).intersect(bedtool, wo=True):
            dist = int(i[2])-int(i[1])
            tad = tuple(i[:3])
            lfc_list.append(float(i[9]))
        val_dict[u] = lfc_list

    val_dict[10] = [float(x) for x in list(zip(*peak_genes))[6]]

    test_df = make_df_from_dict(val_dict)
    fig, ax = plt.subplots()
    sns.barplot(x="labels", y="values", data=test_df, ax=ax)
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax)  
    ax.set_ylim([-4, 4])
    fig.savefig(f'./plots/tad_interactions/chrom_chipseq/{chrom}_diffRNA_atac.png')
    plt.close(fig)


def compartment_barplot(chrom, chrom_cluster_dict, compartment_vector_dict):
    val_dict = {}
    clusts = chrom_cluster_dict[chrom]
    for u in np.unique(clusts):
        inds = np.where(clusts == u)[0]
        lfc_list = compartment_vector_dict[chrom][inds]
        val_dict[u] = list(lfc_list)

    test_df = make_df_from_dict(val_dict)
    fig, ax = plt.subplots()
    sns.barplot(x="labels", y="values", data=test_df, ax=ax)
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax)  
    ax.set_ylim([-.03, .03])
    fig.savefig(f'./plots/tad_interactions/chrom_chipseq/{chrom}_compartment_atac.png')
    plt.close(fig)


## Makes barplots for individual chromosome clusters
def chromosome_specific_barplot(chrom, chrom_cluster_dict, compartment_vector_dict, peak_genes, 
    index_to_tad_merged, bedtool_dict, atac_count,  all_tcon_v_treg):
    val_dict = {}
    clusts = chrom_cluster_dict[chrom]
    for u in np.unique(clusts):
        inds = np.where(clusts == u)[0]
        lfc_list = compartment_vector_dict[chrom][inds]
        val_dict[u] = list(lfc_list)

    test_df = make_df_from_dict(val_dict)

    n  = len(bedtool_dict)*2 + 1 + 1 + 1
    fig, ax = plt.subplots(n//4+1, 4, figsize=(20, 3*n//4))
    ax = np.ravel(ax)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[0])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[0])
    ax[0].set_ylim([-.03, .03])


    bedtool = peak_genes
    val_dict = {}
    clusts = chrom_cluster_dict[chrom]
    for u in np.unique(clusts):
        lfc_list = []
        inds = np.where(clusts == u)[0]
        changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
        for i in pbt.BedTool(changing_tads).intersect(bedtool, wo=True):
            dist = int(i[2])-int(i[1])
            tad = tuple(i[:3])
            lfc_list.append(float(i[11]))
        val_dict[u] = lfc_list

    val_dict[10] = [float(x) for x in list(zip(*peak_genes))[6]]

    test_df = make_df_from_dict(val_dict)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[1])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[1])
    ax[1].set_ylim([0, 5_000])

    ax_counter = 1
    # Normalized by ATAC-seq
    val_dict = {}
    for name in bedtool_dict:
        bedtool = bedtool_dict[name]
        clusts = chrom_cluster_dict[chrom]
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
            
            for i in pbt.BedTool(changing_tads).intersect(bedtool.intersect(all_tcon_v_treg, u=True), c=True):
                tad = tuple(i[:3])
                count = int(i[3])
                if atac_count[tad] > 0:
                    lfc_list.append(count/(atac_count[tad]))
                else:
                    continue
                    lfc_list.append(0)
            val_dict[u] = lfc_list
        test_df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=test_df, ax=ax[2*ax_counter])
        sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[2*ax_counter])

        ax[2*ax_counter].set_title(f'{name}_atac')
        top = np.max([np.mean(x) for x in val_dict.values()])
        ax[2*ax_counter].set_ylim([0, 2*top])
        ax_counter += 1

    ax_counter = 1
    val_dict = {}
    for name in bedtool_dict:
        bedtool = bedtool_dict[name]
        clusts = chrom_cluster_dict[chrom]
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
            
            for i in pbt.BedTool(changing_tads).intersect(bedtool.intersect(all_tcon_v_treg, u=True), c=True):
                dist = (int(i[2])-int(i[1]))/1000                
                tad = tuple(i[:3])
                count = int(i[3])
                count = count/dist
                lfc_list.append(count)
            val_dict[u] = lfc_list
        test_df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=test_df, ax=ax[2*ax_counter+1])
        sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[2*ax_counter+1])

        ax[2*ax_counter+1].set_title(name)
        top = np.max([np.mean(x) for x in val_dict.values()])
        ax[2*ax_counter+1].set_ylim([0, 2*top])
        ax_counter +=1
    return fig


def rna_barplot(clusts, compartment_vector, peak_genes, 
    all_tads, bedtool_dict, atac_count,  all_tcon_v_treg):    
    bedtool = peak_genes
    val_dict = {}
    fig, ax = plt.subplots(2, 2, figsize=(12, 6) )
    ax = np.ravel(ax)
    for u in np.unique(clusts):
        lfc_list = []
        inds = np.where(clusts == u)[0]
        changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
        for i in pbt.BedTool(changing_tads).intersect(bedtool, wo=True):
            dist = int(i[2])-int(i[1])
            tad = tuple(i[:3])
            lfc_list.append(float(i[11]))
        val_dict[u] = lfc_list
    val_dict[10] = [float(x) for x in list(zip(*peak_genes))[8]]
    test_df = make_df_from_dict(val_dict)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[0])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[0])
    ax[0].set_ylim([0, 5_000])
    ax[0].set_title("All peak genes")

    val_dict = {}
    for u in np.unique(clusts):
        lfc_list = []
        inds = np.where(clusts == u)[0]
        changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
        for i in pbt.BedTool(changing_tads).intersect(bedtool, wo=True):
            dist = int(i[2])-int(i[1])
            tad = tuple(i[:3])
            lfc_list.append(float(i[9]))
        val_dict[u] = lfc_list
    val_dict[10] = [float(x) for x in list(zip(*peak_genes))[6]]
    test_df = make_df_from_dict(val_dict)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[1])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[1])
    ax[1].set_title("All peak genes")
    ax[1].set_ylim([-6, 6])



    bedtool = peak_genes.filter(lambda x: float(x[4]) < .05).saveas()
    val_dict = {}
    for u in np.unique(clusts):
        lfc_list = []
        inds = np.where(clusts == u)[0]
        changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
        for i in pbt.BedTool(changing_tads).intersect(bedtool, wo=True):
            dist = int(i[2])-int(i[1])
            tad = tuple(i[:3])
            lfc_list.append(float(i[11]))
        val_dict[u] = lfc_list
    val_dict[10] = [float(x) for x in list(zip(*peak_genes))[8]]
    test_df = make_df_from_dict(val_dict)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[2])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[2])
    ax[2].set_title("Peak genes < .05")
    ax[2].set_ylim([0, 5_000])

    for u in np.unique(clusts):
        lfc_list = []
        inds = np.where(clusts == u)[0]
        changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
        for i in pbt.BedTool(changing_tads).intersect(bedtool, wo=True):
            dist = int(i[2])-int(i[1])
            tad = tuple(i[:3])
            lfc_list.append(float(i[9]))
        val_dict[u] = lfc_list
    val_dict[10] = [float(x) for x in list(zip(*peak_genes))[6]]
    test_df = make_df_from_dict(val_dict)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[3])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[3])
    ax[3].set_title("Peak genes < .05")
    ax[3].set_ylim([-6, 6])
    for a in ax:
        a.set_ylabel("")
        a.set_xlabel("")

    return fig




## Makes barplots for final, all-chromosome clusters
def full_barplot(clusts, compartment_vector, peak_genes, 
    all_tads, bedtool_dict, atac_count,  all_tcon_v_treg):
    val_dict = {}
    for u in np.unique(clusts):
        inds = np.where(clusts == u)[0]
        lfc_list = compartment_vector[inds]
        val_dict[u] = list(lfc_list)
    test_df = make_df_from_dict(val_dict)

    n  = len(bedtool_dict)*2 + 1 + 1 + 1
    fig, ax = plt.subplots(n//4+1, 4, figsize=(20, 3*n//4))
    ax = np.ravel(ax)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[0])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[0])
    ax[0].set_ylim([-.03, .03])
    ax[0].set_title("Compartment scores")


    bedtool = peak_genes
    val_dict = {}
    for u in np.unique(clusts):
        lfc_list = []
        inds = np.where(clusts == u)[0]
        changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
        for i in pbt.BedTool(changing_tads).intersect(bedtool, wo=True):
            dist = int(i[2])-int(i[1])
            tad = tuple(i[:3])
            lfc_list.append(float(i[11]))
        val_dict[u] = lfc_list

    val_dict[10] = [float(x) for x in list(zip(*peak_genes))[6]]

    test_df = make_df_from_dict(val_dict)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[1])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[1])
    ax[1].set_ylim([0, 5_000])
    ax[1].set_title("RPKM scores")

    ax_counter = 2
    # Normalized by ATAC-seq
    val_dict = {}
    for name in bedtool_dict:
        bedtool = bedtool_dict[name]
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
            
            for i in pbt.BedTool(changing_tads).intersect(bedtool.intersect(all_tcon_v_treg, u=True), c=True):
                tad = tuple(i[:3])
                count = int(i[3])
                if atac_count[tad] > 0:
                    lfc_list.append(count/(atac_count[tad]))
                else:
                    continue
                    lfc_list.append(0)
            val_dict[u] = lfc_list
        test_df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=test_df, ax=ax[ax_counter])
        sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[ax_counter])

        ax[ax_counter].set_title(f'{name}_atac')
        top = np.max([np.mean(x) for x in val_dict.values()])
        ax[ax_counter].set_ylim([0, 2*top])
        ax_counter += 1

    val_dict = {}
    for name in bedtool_dict:
        bedtool = bedtool_dict[name]
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
            for i in pbt.BedTool(changing_tads).intersect(bedtool, c=True):
                dist = (int(i[2])-int(i[1]))/1000                
                tad = tuple(i[:3])
                count = int(i[3])
                count = count/dist
                lfc_list.append(count)
            val_dict[u] = lfc_list
        test_df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=test_df, ax=ax[ax_counter])
        sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[ax_counter])

        ax[ax_counter].set_title(name)
        top = np.max([np.mean(x) for x in val_dict.values()])
        ax[ax_counter].set_ylim([0, 2*top])
        ax_counter +=1
    for a in ax:
        a.set_ylabel("")
        a.set_xlabel("")
    plt.tight_layout()        
    return fig

def full_bigwig_barplot_with_bbis(file_dict, compartment_vector_dict, clusts, all_tads, peaks, normalize):
    val_dict = {}

    for u in np.unique(clusts):
        inds = np.where(clusts == u)[0]
        lfc_list = compartment_vector_dict[inds]
        val_dict[u] = list(lfc_list)
    test_df = make_df_from_dict(val_dict)
    n  = len(file_dict) + 1
    fig, ax = plt.subplots(n//4+1, 4, figsize=(15, 3*n//4))
    ax = np.ravel(ax)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[0])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[0])
    ax[0].set_ylim([-.03, .03])
    
    ax_count = 1
    for name in file_dict:
        tcon_ac_bw = file_dict[name]
        val_dict = {}
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]

            chroms, ss, es = [], [],[]                
            for i in peaks.intersect(pbt.BedTool(changing_tads), u=True):
                chrom, s, e = i[:3]
                chroms.append('chr' + chrom)
                ss.append(s)
                es.append(e)
            if len(chroms) == 0:
                val_dict[u] = [0]
                continue
            tcon = tcon_ac_bw.stackup(chroms, ss, es, bins=10000)
            norm = normalize.stackup(chroms, ss, es, bins=10000)
            tcon_ac_scores = np.nanmean(tcon - norm, axis=1)[:, None]
            fcs = np.squeeze((tcon_ac_scores), axis=1)
            val_dict[u] = list(fcs)
        df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=df, ax=ax[ax_count])
        sns.stripplot(x="labels", y="values", data=df, color="0", alpha=.5, ax=ax[ax_count])  
        ax[ax_count].set_title(name)
        ax[ax_count].set_ylim([-5, 5])
        ax_count += 1

    for a in ax:
        a.set_xlabel("")
        a.set_ylabel("")
    plt.tight_layout()
    return fig




def full_bigwig_barplot(file_dict, compartment_vector_dict, clusts, all_tads):
    val_dict = {}

    for u in np.unique(clusts):
        inds = np.where(clusts == u)[0]
        lfc_list = compartment_vector_dict[inds]
        val_dict[u] = list(lfc_list)
    test_df = make_df_from_dict(val_dict)
    n  = len(file_dict) + 1
    fig, ax = plt.subplots(n//2+1, 4, figsize=(15, 3*n//2))
    ax = np.ravel(ax)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[0])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[0])
    ax[0].set_ylim([-.03, .03])
    

    ax_count = 1
    for name in file_dict:
        file_treg, file_tcon = file_dict[name]
        tcon_ac_bw = bbi.open(f'/Users/gabrieladolsten/Downloads/{file_tcon}')
        val_dict = {}
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]

            chroms, ss, es = [], [],[]                
            for i in pbt.BedTool(changing_tads):
                chrom, s, e = i[:3]
                chroms.append('chr' + chrom)
                ss.append(s)
                es.append(e)
            tcon_ac_scores = np.nanmean(tcon_ac_bw.stackup(chroms, ss, es, bins=10000), axis=1)[:, None]
            fcs = np.squeeze((tcon_ac_scores), axis=1)
            val_dict[u] = list(fcs)
        df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=df, ax=ax[2*ax_count])
        sns.stripplot(x="labels", y="values", data=df, color="0", alpha=.5, ax=ax[2*ax_count])  
        ax[2*ax_count].set_title(name + " TCON")
        ax[2*ax_count].set_ylim([0, 2])
        ax_count += 1

    ax_count = 1
    for name in file_dict:
        file_treg, file_tcon = file_dict[name]
        treg_ac_bw = bbi.open(f'/Users/gabrieladolsten/Downloads/{file_treg}')
        val_dict = {}
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(all_tads)[inds, :]]
            chroms, ss ,es, = [], [],[]                
            for i in pbt.BedTool(changing_tads):
                chrom, s, e = i[:3]
                chroms.append('chr' + chrom)
                ss.append(s)
                es.append(e)
            treg_ac_scores = np.nanmean(treg_ac_bw.stackup(chroms, ss, es, bins=10000), axis=1)[:, None]
            fcs = np.squeeze((treg_ac_scores), axis=1)
            val_dict[u] = list(fcs)
        df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=df, ax=ax[2*ax_count+1])
        sns.stripplot(x="labels", y="values", data=df, color="0", alpha=.5, ax=ax[2*ax_count+1])  
        ax[2*ax_count+1].set_title(name + " TREG")
        ax[2*ax_count+1].set_ylim([0, 2])
        ax_count += 1
    for a in ax:
        a.set_xlabel("")
        a.set_ylabel("")
    plt.tight_layout()
    return fig


def bigwig_barplot(file_dict, compartment_vector_dict, chrom, chrom_cluster_dict, index_to_tad_merged):
    val_dict = {}
    clusts = chrom_cluster_dict[chrom]
    for u in np.unique(clusts):
        inds = np.where(clusts == u)[0]
        lfc_list = compartment_vector_dict[chrom][inds]
        val_dict[u] = list(lfc_list)

    test_df = make_df_from_dict(val_dict)
    n  = len(file_dict) + 1
    fig, ax = plt.subplots(n//2+1, 4, figsize=(15, 3*n//2))
    ax = np.ravel(ax)
    sns.barplot(x="labels", y="values", data=test_df, ax=ax[0])
    sns.stripplot(x="labels", y="values", data=test_df, color="0", alpha=.1, ax=ax[0])
    ax[0].set_ylim([-.03, .03])
    

    ax_count = 1
    for name in file_dict:
        file_treg, file_tcon = file_dict[name]
        tcon_ac_bw = bbi.open(f'/Users/gabrieladolsten/Downloads/{file_tcon}')
        val_dict = {}
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
            chroms, ss ,es, = [], [],[]                
            for i in pbt.BedTool(changing_tads):
                chrom, s, e = i[:3]
                chroms.append('chr' + chrom)
                ss.append(s)
                es.append(e)
            tcon_ac_scores = np.nanmean(tcon_ac_bw.stackup(chroms, ss, es, bins=10000), axis=1)[:, None]
            fcs = np.squeeze((tcon_ac_scores), axis=1)
            val_dict[u] = list(fcs)
        df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=df, ax=ax[2*ax_count])
        sns.stripplot(x="labels", y="values", data=df, color="0", alpha=.5, ax=ax[2*ax_count])  
        ax[2*ax_count].set_title(name + " TCON")
        ax[2*ax_count].set_ylim([0, 3])
        ax_count += 1

    ax_count = 1
    for name in file_dict:
        file_treg, file_tcon = file_dict[name]
        treg_ac_bw = bbi.open(f'/Users/gabrieladolsten/Downloads/{file_treg}')
        val_dict = {}
        for u in np.unique(clusts):
            lfc_list = []
            inds = np.where(clusts == u)[0]
            changing_tads = [list(x) for x in np.asarray(index_to_tad_merged[chrom])[inds, :]]
            chroms, ss ,es, = [], [],[]                
            for i in pbt.BedTool(changing_tads):
                chrom, s, e = i[:3]
                chroms.append('chr' + chrom)
                ss.append(s)
                es.append(e)
            treg_ac_scores = np.nanmean(treg_ac_bw.stackup(chroms, ss, es, bins=10000), axis=1)[:, None]
            fcs = np.squeeze((treg_ac_scores), axis=1)
            val_dict[u] = list(fcs)
        df = make_df_from_dict(val_dict)
        sns.barplot(x="labels", y="values", data=df, ax=ax[2*ax_count+1])
        sns.stripplot(x="labels", y="values", data=df, color="0", alpha=.5, ax=ax[2*ax_count+1])  
        ax[2*ax_count+1].set_title(name + " TREG")
        ax[2*ax_count+1].set_ylim([0, 3])
        ax_count += 1
    for a in ax:
        a.set_xlabel("")
        a.set_ylabel("")
    plt.tight_layout()
    return fig



def make_all_compartment_vectors(parsed_chroms, tad_to_index_merged, tn_compartment, treg_compartment):
    compartment_vector_dict = {}
    for chrom in parsed_chroms:
        compartment_vector = make_compartment_vector(chrom, tad_to_index_merged, tn_compartment)
        compartment_vector_dict[chrom] = compartment_vector

    compartment_interchrom_vector_dict = {}
    for chrom in parsed_chroms:
        compartment_vector = []
        for chrom2 in parsed_chroms:
            if chrom == chrom2:
                continue
            else:
                compartment_vector += list(make_compartment_vector(chrom2, tad_to_index_merged, tn_compartment))
        compartment_interchrom_vector_dict[chrom] = np.asarray(compartment_vector)

        
    treg_compartment_vector_dict = {}
    for chrom in parsed_chroms:
        compartment_vector = make_compartment_vector(chrom, tad_to_index_merged, treg_compartment)
        treg_compartment_vector_dict[chrom] = compartment_vector

    treg_compartment_interchrom_vector_dict = {}
    for chrom in parsed_chroms:
        compartment_vector = []
        for chrom2 in parsed_chroms:
            if chrom == chrom2:
                continue
            else:
                compartment_vector += list(make_compartment_vector(chrom2, tad_to_index_merged, treg_compartment))
        treg_compartment_interchrom_vector_dict[chrom] = np.asarray(compartment_vector)
        
        
    size_vector_dict = {}
    for chrom in parsed_chroms:
        sizes = [int(x[2]) - int(x[1]) for x in tad_to_index_merged[chrom]]
        size_vector_dict[chrom] = sizes
        
    size_interchrom_vector_dict = {}
    for chrom in parsed_chroms:
        all_sizes = [] 
        for chrom2 in parsed_chroms:
            if chrom == chrom2:
                continue
            else:
                sizes = [int(x[2]) - int(x[1]) for x in tad_to_index_merged[chrom2]]
            all_sizes += sizes
        size_interchrom_vector_dict[chrom] = np.asarray(all_sizes)    
    return compartment_vector_dict, compartment_interchrom_vector_dict, treg_compartment_vector_dict, treg_compartment_interchrom_vector_dict, size_vector_dict, size_interchrom_vector_dict
