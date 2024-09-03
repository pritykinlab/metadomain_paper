import numpy as np
import matplotlib.pyplot as plt
import pyBigWig
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
import seaborn as sns


def make_delta_chip_and_looping_df(valdf, has_treg_specific, has_tcon_specific):
    # For each loop which is gained, check which cluster
    df = pd.DataFrame()
    for i in zip(*np.where(has_treg_specific & has_tcon_specific)):
        treg_specific_loops = np.where(treg_specific[i])[0]
        tcon_specific_loops = np.where(tcon_specific[i])[0]
        treg_subdf = valdf[valdf['Inds'].isin(treg_specific_loops)]
        tcon_subdf = valdf[valdf['Inds'].isin(tcon_specific_loops)]

        ind_subdf = valdf[valdf['Inds'].isin(i)]
        h3k27ac_delta = np.nanmean(ind_subdf['h3k27ac_treg']-ind_subdf['h3k27ac_tcon'])
        h3k27me3_delta = np.nanmean(ind_subdf['h3k27me3_treg']-ind_subdf['h3k27me3_tcon'])

        treg_ind_tcon_27ac = np.nanmean(treg_subdf['h3k27ac_tcon'])
        tcon_ind_tcon_27ac = np.nanmean(tcon_subdf['h3k27ac_tcon'])

        treg_ind_tcon_27me3 = np.nanmean(treg_subdf['h3k27me3_tcon'])
        tcon_ind_tcon_27me3 = np.nanmean(tcon_subdf['h3k27me3_tcon'])

        row = pd.DataFrame([[h3k27ac_delta, h3k27me3_delta, treg_ind_tcon_27ac, 
         tcon_ind_tcon_27ac, treg_ind_tcon_27me3, tcon_ind_tcon_27me3]])
        row.index = i
        df = pd.concat([df, row], axis=0)
    df.columns = ['∆ H3K27ac', '∆ H3K27me3', 'H3K27ac on Treg inds', 
         'H3K27ac on Tcon inds', 'H3K27me3 on Treg inds', 'H3K27me3 on Tcon inds']
    df.index=np.where(has_treg_specific & has_tcon_specific)[0]
    return df


def plot_change_in_loops_due_to_chip(df_treg_all, df_tcon_all, delta_col = 'H3K27ac', delta_co = .5, chip_cluster_map=None):
	sns.set_style("whitegrid")
	indsoi = df_treg_all[f'∆ {delta_col}'] > delta_co
	subdf = df_treg_all[indsoi].sum(axis=0)-df_tcon_all[indsoi].sum(axis=0)[2:]
	subdf = subdf[:-2]
	subdf = subdf[['H3K9me3', f'H3K27ac', 'H3K27me3']]
	subdf = pd.DataFrame(subdf)
	subdf['Condition'] = f'Treg {delta_col} > Tcon'

	indsoi = df_treg_all[f'∆ {delta_col}'] < -delta_co
	subdf2 = df_treg_all[indsoi].sum(axis=0)-df_tcon_all[indsoi].sum(axis=0)[2:]
	subdf2 = subdf2[:-2]
	subdf2 = subdf2[['H3K9me3', f'H3K27ac', 'H3K27me3']]
	subdf2 = pd.DataFrame(subdf2)
	subdf2['Condition'] = f'Tcon {delta_col} > Treg'

	tmpdf = pd.concat([subdf, subdf2], axis=0)
	tmpdf = tmpdf.reset_index()
	tmpdf.columns = ['ChIP Cluster', 'Change', 'Condition']
	tmpdf = tmpdf.melt(["Condition", "ChIP Cluster"], 'Change')

	fig, ax = plt.subplots(figsize=(4, 4))
	sns.barplot(data=tmpdf, x='ChIP Cluster', y='value', 
	            hue='Condition', palette=['red', 'blue'],)
	#            hue_colors=['red', "blue"])
	leg = ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
	leg.set_title("Bins where:")
	ax.set_ylabel("# Megaloops")
	ax.set_title("Change in number of \nmegaloops with \nChIP Cluster (Treg - Tcon)")
	ax.set_xlabel("ChIP Cluster")
	plt.xticks(rotation=20)
	ax.set_yticks([-50, -25, 0, 25, 50])
	ax.set_xticks([0, 1, 2])
	ax.set_ylim([-50, 50])
	plt.xticks(y=-.1)
	newax = ax.inset_axes(transform=ax.transAxes,
		   bounds = (0, -.1, 1, .1))
	newax.matshow(np.arange(3)[None, :], cmap=chip_cluster_map, aspect='auto')
	newax.set_xticks([])
	newax.set_yticks([])



	return fig, ax, subdf, subdf2

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

import seaborn as sns
import scipy
import scipy.ndimage
def plot_valdf_chipseq_clusters(comp_vec, valdf, linkage, compartment_label='Compartment (Tcon)', title='ChIP-seq over\nmegaloop bins',
				chip_cluster_map = None, keys=None, col2label={}, comp_cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)):
	sns.set(font_scale=1.5)
	sns.set_style('whitegrid')
	comp_scores = comp_vec[valdf['Inds']]
	subdf = valdf[keys].copy()
	clust = valdf['ChIP_Clust'].values
	tcon_col_ind = arr(['treg' in x for x in subdf.columns])
	subdf.columns = [col2label[x] for x in subdf.columns]
	cmap = LinearSegmentedColormap.from_list(
	    name='test', 
	    colors=['white', 'brown']
	)

	# blue = 0x44, 0x6D, 0x95
	# green = 0x13, 0x5F, 0x0F
	# amber = 0xD6, 0x8B, 0x1A

	f = sns.clustermap(subdf, cmap=cmap, vmin=1, vmax=3, 
	            method = 'average', metric='cosine', figsize=(4, 5.5), row_linkage = linkage,
	                cbar_pos=(1.0, .63, 0.05, 0.18), rasterized=True, xticklabels=True)
	ax = f.ax_heatmap
	f.ax_col_dendrogram.remove()
	ax.set_yticks([])
	ax.set_title(title, y=1.15)
	ax.set_ylabel("Megaloop bins")
	ax.yaxis.set_label_coords(-.6, .5)
	plt.sca(ax)
	plt.xticks(fontsize=12)

	_ = f.dendrogram_row.reordered_ind

	d = .2
	newax = ax.inset_axes(transform=ax.transAxes,
	       bounds = (1+d*1+d/3, 0, d*2/3, 1))
	newax.matshow(comp_scores[_][:, None], cmap=comp_cmap, 
	              vmin=-1.5, vmax=1.5, aspect='auto')
	newax.set_xticks([])
	newax.set_yticks([])
	newax.set_zorder(0)
	newax.set_title(compartment_label, rotation=40, ha='left', x=.3, fontsize=12)

	newax = ax.inset_axes(transform=ax.transAxes,
	       bounds = (1+d*0+d/3, 0, d*2/3, 1))
	newax.matshow(clust[_][:, None], cmap=chip_cluster_map, aspect='auto',)
	newax.set_xticks([])
	newax.set_yticks([])
	newax.set_zorder(0)
	newax.set_title("ChIP Clustering", rotation=40, ha='left', x=.1, fontsize=12, )

	boo = (valdf['n_MegaLoops'] > 1).values[_]
	#boo = scipy.ndimage.gaussian_filter1d(boo[_].astype(float), sigma=1) > .1
	newax = ax.inset_axes(transform=ax.transAxes,
		   bounds = (1+d*2+d/3, 0, d*2/3, 1))
	newax.matshow(boo[:, None], cmap=cm.bwr,
				  vmin=-2, vmax=2, aspect='auto')
	newax.set_xticks([])
	newax.set_yticks([])
	newax.set_zorder(0)
	newax.set_title("Megaloop", rotation=40, ha='left', x=.3, fontsize=12); assert len(boo) == len(comp_scores); assert len(boo) == len(subdf);

	_col = f.dendrogram_col.reordered_ind; print(len(boo), len(comp_scores), len(tcon_col_ind), len(subdf))
	newax = ax.inset_axes(transform=ax.transAxes,
	       bounds = (0, -.04, 1, .04))
	newax.matshow(tcon_col_ind[_col][None, :], cmap=cm.bwr, aspect='auto',); print(len(tcon_col_ind))
	newax.set_xticks([]); newax.set_yticks([])
	return f, _



def plot_enrichment_df(df, int_to_chip_clust, chip_cluster_map):
	fig, axs = plt.subplots(1, 2, figsize=(10, 4))
	ax = axs[0]
	axs[0].remove()

	cmap = LinearSegmentedColormap.from_list(
	    name='test', 
	    colors=['lightgray', 'white', 'brown']
	)

	jointloop_mat = df.values; n = len(jointloop_mat)
	fraction_loops1 = jointloop_mat.sum(axis=0)
	fraction_loops2 = jointloop_mat.sum(axis=1)
	fraction_loops2 = fraction_loops2/fraction_loops2.sum()
	expected_mat = np.outer(fraction_loops1, fraction_loops2)

	enrichment_loop_mat = pd.DataFrame(np.log2(jointloop_mat/expected_mat))

	enrichment_loop_mat.columns = pd.Series(enrichment_loop_mat.columns).apply(lambda x: int_to_chip_clust[x])
	enrichment_loop_mat.index = df.columns

	ax = axs[1]
	sns.heatmap(enrichment_loop_mat, cmap=cmap, vmin=-1, vmax=1, annot=True, ax=ax, square=True)
	ax.set_title("Enrichment for megaloops \nbetween ChIP clusters", y=1.1)
	plt.yticks(x=-.1)
	plt.xticks(y=-.1)
	newax = ax.inset_axes(transform=ax.transAxes,
	       bounds = (-.1, 0, .1, 1))
	newax.matshow(np.arange(n)[:, None], cmap=chip_cluster_map, aspect='auto')
	newax.set_xticks([])
	newax.set_yticks([])

	newax = ax.inset_axes(transform=ax.transAxes,
	       bounds = (0, -.1, 1, .1))
	newax.matshow(np.arange(n)[None, :], cmap=chip_cluster_map, aspect='auto')
	newax.set_xticks([])
	newax.set_yticks([])
	newax.set_zorder(0)
	plt.xticks(rotation=20)
	plt.yticks(rotation=0)
	ax.set_xlabel('Cluster of megaloop partner')
	ax.set_ylabel('Cluster of \nmegaloop bin', rotation=90)
	ax.yaxis.set_label_coords(-.65, .5)
	return fig, ax, enrichment_loop_mat

from aux_functions import *
import pybedtools as pbt
import sklearn
from sklearn.decomposition import PCA
from scipy.stats import zscore
def stack_chip_over_megaloops(full_cluster_df, bwdict, all_ind_to_region, keys=['h3k27ac_tcon', 'h3k27me3_tcon', 'h3k27ac_treg', 'h3k27me3_treg'],
                                valmin=0, valmax=np.inf, cluster_method='ward', cluster_metric='euclidean', n_clusters=3):
        valdf = pd.DataFrame()
        for key in keys:
            inds = full_cluster_df['Inds']
            regs = pbt.BedTool([all_ind_to_region[x] for x in inds])
            chrom, s, e = unzreg(regs)
            chrom = add_chr_to_list(chrom)
            v = np.nanmean(bwdict[key].stackup(chrom, s, e, bins=200), axis=1)
            valdf[key] = v

        valdf['Inds'] = full_cluster_df['Inds'].values
        valdf['goodClust'] = full_cluster_df['goodClust'].values
        valdf['Clusts'] = full_cluster_df['Clust'].values
        valmat = np.clip(valmin, valmax, valdf[keys].values)
        valmat = np.log(valmat+.5)
        #valmat = zscore(valmat)
        pca = PCA(n_components=2)
        pca = pca.fit(valmat)
        components = pca.transform(valmat); print(components.shape)
        _, clust, linkage = make_order_and_cluster_custom(components, n_clusters=n_clusters,
                                                 method = cluster_method, metric=cluster_metric)
        valdf['ChIP_Clust'] = rename_clusts_by_order(clust, _)
        valdf.index = valdf['Inds']
        return valdf, linkage


import scipy
def make_order_and_cluster_custom(matrix, method='average', metric='cosine', n_clusters=2):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric=metric)
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)
        order = dendro['leaves']
        poop = scipy.cluster.hierarchy.fcluster(linkage, t=n_clusters, criterion='maxclust')-1
    else:
        order = np.arange(len(matrix))
    return order, poop, linkage

def make_int_to_chip_clust(valdf):
        grouped_df = valdf[['h3k27ac_tcon', 'h3k27me3_tcon', 'h3k9me3_tcon', 'h3k27ac_treg', 'h3k27me3_treg', 'ChIP_Clust']].groupby('ChIP_Clust').mean()
        h3k27ac_cluster = grouped_df['h3k27ac_tcon'].argmax()
        h3k27me3_cluster = grouped_df['h3k27me3_tcon'].argmax()
        h3k9me3_cluster = grouped_df['h3k9me3_tcon'].argmax()
        misc_cluster = list(set([0, 1, 2, 3])-set([h3k27ac_cluster, h3k27me3_cluster, h3k9me3_cluster]))[0]
        int_to_chip_clust = {}
        int_to_chip_clust[h3k9me3_cluster] = 'H3K9me3'
        int_to_chip_clust[h3k27me3_cluster] = 'H3K27me3'
        int_to_chip_clust[h3k27ac_cluster] = 'H3K27ac'
        int_to_chip_clust[misc_cluster] = 'Misc.'
        #assert len(int_to_chip_clust) == 4
        return int_to_chip_clust


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


def get_bedtool_lengths(bedtool):
    es = get_col(bedtool, 2).astype(int) 
    ss = get_col(bedtool, 1).astype(int)
    deltas = np.abs(es-ss)
    return deltas

def get_basemeans(valdf, bedtool_oi, chip_clust_to_int, my_treg_comp, all_ind_to_region):
    connection_indset = set(valdf[valdf['ChIP_Clust'] == chip_clust_to_int['H3K27ac']]['Inds'].values)
    connection_regs = pbt.BedTool([all_ind_to_region[x] for x in connection_indset])
    genesoi = bedtool_oi.intersect(connection_regs, u=True)
    cnxn_basemeans = get_col(genesoi, -2).astype(float)
    cnxn_lens = get_bedtool_lengths(genesoi)
    cnxn_h3k27ac = cnxn_basemeans/cnxn_lens

    connection_indset = set(valdf[valdf['ChIP_Clust'] == chip_clust_to_int['H3K27me3']]['Inds'].values)
    connection_regs = pbt.BedTool([all_ind_to_region[x] for x in connection_indset])
    genesoi = bedtool_oi.intersect(connection_regs, u=True)
    cnxn_basemeans = get_col(genesoi, -2).astype(float)
    cnxn_lens = get_bedtool_lengths(genesoi)
    cnxn_h3k27me3 = cnxn_basemeans/cnxn_lens

    connection_indset = set(valdf[valdf['ChIP_Clust'] == chip_clust_to_int['Misc.']]['Inds'].values)
    connection_regs = pbt.BedTool([all_ind_to_region[x] for x in connection_indset])
    genesoi = bedtool_oi.intersect(connection_regs, u=True)
    cnxn_basemeans = get_col(genesoi, -2).astype(float)
    cnxn_lens = get_bedtool_lengths(genesoi)
    cnxn_misc = cnxn_basemeans/cnxn_lens


    aset = set(np.where(my_treg_comp>0)[0])
    acomp_regs = pbt.BedTool([all_ind_to_region[x] for x in aset])
    genesoi = bedtool_oi.intersect(acomp_regs, u=True)
    acomp_basemeans = get_col(genesoi, -2).astype(float)
    acomp_lengths = get_bedtool_lengths(genesoi)
    cnxn_acomp = (acomp_basemeans/acomp_lengths)
    d = {
         'H3K27me3 Megaloops' : list(cnxn_h3k27me3),
         'Misc. Megaloops' : list(cnxn_misc),
         'A compartment' : list(cnxn_acomp), 
         'H3K27ac Megaloops' : list(cnxn_h3k27ac),
    }
    return d



def get_rpkm(indsoi, all_ind_to_region, bedtool_oi, ):
    connection_regs = pbt.BedTool([all_ind_to_region[x] for x in indsoi])
    genesoi = bedtool_oi.intersect(connection_regs, u=True)
    cnxn_basemeans = get_col(genesoi, -2).astype(float)
    cnxn_lens = get_bedtool_lengths(genesoi)
    cnxn_rpkm = cnxn_basemeans/cnxn_lens*1_000
    return cnxn_rpkm, genesoi



import networkx as nx
indco = 10
def get_n_of_distinct_regions(inds, indco):
    m = inds.values
    diff = np.abs(np.subtract.outer(m, m)) < indco
    G = nx.from_numpy_matrix(diff)
    conn_comp = nx.connected_components(G)
    l = list(conn_comp)
    
    clusts = np.zeros_like(m)
    for c, i in enumerate(l):
        i = list(i)
        clusts[i]=c
    assert len(l) == len(np.unique(clusts))
    return len(l), clusts

def make_o_df_clust(final_intra_focal_contacts_treg, final_intra_focal_contacts_tcon, all_intra_connections, s, e, n_clusts=7):
    connections = pd.DataFrame(all_intra_connections[s:e, s:e])
    connections.index = np.arange(s, e)
    connections.columns = np.arange(s, e)
    indsoi = connections.sum(axis=1)>2
    subcon = connections.loc[indsoi, :].loc[:, indsoi]
    o, clust = make_order2_and_cluster(subcon, n_clusters=n_clusts)
    
    treg_connections = pd.DataFrame(final_intra_focal_contacts_treg[s:e, s:e])
    treg_connections.index = np.arange(s, e)
    treg_connections.columns = np.arange(s, e)
    treg_subcon = treg_connections.loc[indsoi, :].loc[:, indsoi]

    tcon_connections = pd.DataFrame(final_intra_focal_contacts_tcon[s:e, s:e])
    tcon_connections.index = np.arange(s, e)
    tcon_connections.columns = np.arange(s, e)
    tcon_subcon = tcon_connections.loc[indsoi, :].loc[:, indsoi]
    return subcon, o, clust, treg_subcon, tcon_subcon

def extract_networks(df, clusters, indco=10):
    # assumes mat and clusters ordered the same way
    goodclusts = []
    us = np.unique(clusters)
    frac_different_list = []
    n_of_distinct_regions_list = []
    connection_freq_list = []
    for u in us:
        submat = df.iloc[clusters==u, :].iloc[:, clusters==u]
        inds = submat.index
        n_of_distinct_regions, _ = get_n_of_distinct_regions(inds, indco)
        frac_different = n_of_distinct_regions/len(submat)
        connection_freq = submat.values.mean()
        if (frac_different > .4) and (n_of_distinct_regions > 2) and (connection_freq > .3):
            goodclusts.append(u)
            
        frac_different_list.append(frac_different)
        n_of_distinct_regions_list.append(n_of_distinct_regions)
        connection_freq_list.append(connection_freq)
    return goodclusts, (frac_different_list, n_of_distinct_regions_list, connection_freq_list)





from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from tad_functions import make_order2_and_cluster
def make_cluster_df_and_full_cluster_df(final_intra_focal_contacts_treg, final_intra_focal_contacts_tcon, comp_vec, chrom_to_start, chrom_to_end, all_intra_connections, n_clusts=7):
    labeldict = {
        -1:'Tcon',
         1:'Treg',
         0:'Tcon and Treg',
    }
    sns.set(font_scale=1.5)
    sns.set_style('whitegrid')
    n_clusts = 7
    comp = comp_vec.copy()
    comp[np.isnan(comp)]=0
    inds_to_keep = []
    clusts_of_inds = []
    full_cluster_df = pd.DataFrame(columns=["Inds", "Clust", 'goodClust', "Chrom"])

    cluster_df = pd.DataFrame()
    for chrom in chrom_to_start:
        s, e = chrom_to_start[chrom], chrom_to_end[chrom]

        subcon, o, clust, treg_subcon, tcon_subcon = make_o_df_clust(final_intra_focal_contacts_treg, final_intra_focal_contacts_tcon, all_intra_connections, s, e, n_clusts=n_clusts)
        clust = rename_clusts_by_order(clust, o)
        osubmat = subcon.iloc[:, o].iloc[o, :]; 
        inds = osubmat.index
        clustsoi, diagnostic = extract_networks(osubmat, clust[o])


        inds_with_good_cluster = np.isin(clust[o], clustsoi)
        goodinds = inds[inds_with_good_cluster]
        inds_to_keep += list(goodinds)
        running_clusts = list(arr(clust[o]) + np.max(clusts_of_inds + [0]))

        clusts_of_inds += list(arr(running_clusts)[inds_with_good_cluster])

        tmpdf = pd.DataFrame()
        tmpdf['Inds'] = inds
        tmpdf['Clust'] = running_clusts
        tmpdf['goodClust'] = inds_with_good_cluster
        tmpdf['Chrom'] = [chrom]*len(inds)
        
        full_cluster_df = full_cluster_df.append(tmpdf)
        
        assert (((treg_subcon + tcon_subcon)>0) == subcon).all(axis=None)
        subcon = subcon.astype(float)
        subcon[treg_subcon&tcon_subcon] = 0
        subcon[treg_subcon&(~tcon_subcon)] = 1
        subcon[(~treg_subcon)&tcon_subcon] = -1
        subcon[(~treg_subcon)&(~tcon_subcon)] = np.nan    

        osubmat = subcon.iloc[:, o].iloc[o, :].values.astype(float)

    cluster_df['Inds'] = inds_to_keep
    cluster_df['Clusts'] = clusts_of_inds
    full_cluster_df['Inds'] = full_cluster_df['Inds'].astype(int)
    full_cluster_df['Clust'] = full_cluster_df['Clust'].astype(int)

    return cluster_df, full_cluster_df

def plot_cluster_df_and_full_cluster_df(final_intra_focal_contacts_treg, final_intra_focal_contacts_tcon, 
    comp_vec, chrom_to_start, chrom_to_end, all_intra_connections, valdf, n_clusts=7, comp_cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True),
    chip_cluster_map=None):
    labeldict = {
        -1:'Tcon',
         1:'Treg',
         0:'Tcon and Treg',
    }
    sns.set(font_scale=1.5)
    sns.set_style('whitegrid')
    n_clusts = 7
    comp = comp_vec.copy()
    comp[np.isnan(comp)]=0
    inds_to_keep = []
    clusts_of_inds = []
    full_cluster_df = pd.DataFrame(columns=["Inds", "Clust", 'goodClust', "Chrom"])

    cluster_df = pd.DataFrame()
    for chrom in chrom_to_start:
        s, e = chrom_to_start[chrom], chrom_to_end[chrom]

        subcon, o, clust, treg_subcon, tcon_subcon = make_o_df_clust(final_intra_focal_contacts_treg, final_intra_focal_contacts_tcon, all_intra_connections, s, e, n_clusts=n_clusts)
        clust = rename_clusts_by_order(clust, o)
        osubmat = subcon.iloc[:, o].iloc[o, :]; 
        inds = osubmat.index
        clustsoi, diagnostic = extract_networks(osubmat, clust[o])


        inds_with_good_cluster = np.isin(clust[o], clustsoi)
        goodinds = inds[inds_with_good_cluster]
        inds_to_keep += list(goodinds)
        running_clusts = list(arr(clust[o]) + np.max(clusts_of_inds + [0]))

        clusts_of_inds += list(arr(running_clusts)[inds_with_good_cluster])

        tmpdf = pd.DataFrame()
        tmpdf['Inds'] = inds
        tmpdf['Clust'] = running_clusts
        tmpdf['goodClust'] = inds_with_good_cluster
        tmpdf['Chrom'] = [chrom]*len(inds)
        
        full_cluster_df = full_cluster_df.append(tmpdf)
        

        fig, axs = plt.subplots(1, 2, figsize=(8, 6))
        ax = axs[0]

    #     ax.text(0, 0, s=f'n={len(osubmat)}', transform=ax.transAxes)
        
        

        chip_clusters = valdf.loc[inds]['ChIP_Clust'].values    
        newax = ax.inset_axes(transform=ax.transAxes,
               bounds = (-.12, 0, .05, 1))
        newax.matshow(chip_clusters[:, None], cmap=chip_cluster_map, aspect='auto')
        newax.set_title("ChIP Clustering", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)
        newax.set_xticks([])
        newax.set_yticks([])

        newax = ax.inset_axes(transform=ax.transAxes,
               bounds = (-.24, 0, .05, 1))
        
        

        newax.matshow(comp[inds][:, None], cmap=comp_cmap, vmin=-1, vmax=1, aspect='auto')
        newax.set_title("A Compartment", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)
        newax.set_xticks([])
        newax.set_yticks([])
        
        newax = ax.inset_axes(transform=ax.transAxes,
               bounds = (-.36, 0, .05, 1))
        newax.matshow(inds_with_good_cluster[:, None], cmap=cm.gist_heat_r, vmin=0, vmax=2, aspect='auto')
        newax.set_title("Hub-like", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)
        newax.set_xticks([])
        newax.set_yticks([])

        
        cmap = LinearSegmentedColormap.from_list(
            name='test', 
            colors=['blue', 'white', 'lightgray','white','red',]
        )    
        ax = axs[0]
        assert (((treg_subcon + tcon_subcon)>0) == subcon).all(axis=None)
        subcon = subcon.astype(float)
        subcon[treg_subcon&tcon_subcon] = 0
        subcon[treg_subcon&(~tcon_subcon)] = 1
        subcon[(~treg_subcon)&tcon_subcon] = -1
        subcon[(~treg_subcon)&(~tcon_subcon)] = np.nan    

        osubmat = subcon.iloc[:, o].iloc[o, :].values.astype(float)
        ax.matshow(osubmat, cmap=cmap, vmin=-1, vmax=1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"Bins with megaloops")
        ax.set_xlabel(f'chr{chrom}: n={len(osubmat)}')
        assert len(chip_clusters) == len(osubmat)
        plt.tight_layout()
        us = [0, 1, -1]
        ax.legend([mpatches.Patch(color=cmap((b+1)/2)) for b in us],
               [labeldict[b] for b in us], bbox_to_anchor=(1, 1), loc='upper left', fontsize=12)
        axs[1].remove()
        fig.savefig(f'./plots/figure2/intra_hubs_chr{chrom}.svg')
    cluster_df['Inds'] = inds_to_keep
    cluster_df['Clusts'] = clusts_of_inds
    full_cluster_df['Inds'] = full_cluster_df['Inds'].astype(int)
    full_cluster_df['Clust'] = full_cluster_df['Clust'].astype(int)

    return cluster_df, full_cluster_df




def plot_cluster_df_and_full_cluster_df_V2_NO_VALDF(final_intra_focal_contacts_treg, final_intra_focal_contacts_tcon, 
    comp_vec, chrom_to_start, chrom_to_end, all_intra_connections, n_clusts=7, comp_cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True),
    chip_cluster_map=None, ind_to_gene={}):
    labeldict = {
        -1:'Tcon',
         1:'Treg',
         0:'Tcon and Treg',
    }
    sns.set(font_scale=1.5)
    sns.set_style('ticks')
    n_clusts = 7
    comp = comp_vec.copy()
    comp[np.isnan(comp)]=0
    inds_to_keep = []
    clusts_of_inds = []
    full_cluster_df = pd.DataFrame(columns=["Inds", "Clust", 'goodClust', "Chrom"])

    cluster_df = pd.DataFrame()
    for chrom in chrom_to_start:
        s, e = chrom_to_start[chrom], chrom_to_end[chrom]

        subcon, o, clust, treg_subcon, tcon_subcon = make_o_df_clust(final_intra_focal_contacts_treg, final_intra_focal_contacts_tcon, all_intra_connections, s, e, n_clusts=n_clusts)
        clust = rename_clusts_by_order(clust, o)
        osubmat = subcon.iloc[:, o].iloc[o, :]; 
        inds = osubmat.index
        clustsoi, diagnostic = extract_networks(osubmat, clust[o])


        inds_with_good_cluster = np.isin(clust[o], clustsoi)
        goodinds = inds[inds_with_good_cluster]
        inds_to_keep += list(goodinds)
        running_clusts = list(arr(clust[o]) + np.max(clusts_of_inds + [0]))

        clusts_of_inds += list(arr(running_clusts)[inds_with_good_cluster])

        tmpdf = pd.DataFrame()
        tmpdf['Inds'] = inds
        tmpdf['Clust'] = running_clusts
        tmpdf['goodClust'] = inds_with_good_cluster
        tmpdf['Chrom'] = [chrom]*len(inds)
        
        full_cluster_df = full_cluster_df.append(tmpdf)
        # print(full_cluster_df)
        

        fig, axs = plt.subplots(1, 2, figsize=(8, 6))
        ax = axs[0]
        newax = ax.inset_axes(transform=ax.transAxes,
               bounds = (-.24 + .12, 0, .05, 1))
        newax.matshow(comp[inds][:, None], cmap=comp_cmap, vmin=-1, vmax=1, aspect='auto')
        newax.set_title("Compartment", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)
        newax.set_xticks([])
        newax.set_yticks([])
        
        newax = ax.inset_axes(transform=ax.transAxes,
               bounds = (-.36 + .12, 0, .05, 1))
        newax.matshow(inds_with_good_cluster[:, None], cmap=cm.gist_heat_r, vmin=0, vmax=2, aspect='auto')
        newax.set_title("Hub-like", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)
        newax.set_xticks([])
        newax.set_yticks([])

        
        cmap = LinearSegmentedColormap.from_list(
            name='test', 
            colors=['blue', 'white', 'lightgray','white','red',]
        )    
        ax = axs[0]
        assert (((treg_subcon + tcon_subcon)>0) == subcon).all(axis=None)
        subcon = subcon.astype(float)
        subcon[treg_subcon&tcon_subcon] = 0
        subcon[treg_subcon&(~tcon_subcon)] = 1
        subcon[(~treg_subcon)&tcon_subcon] = -1
        subcon[(~treg_subcon)&(~tcon_subcon)] = np.nan    

        # print(subcon)
        subcon = subcon.iloc[:, o].iloc[o, :]
        genenames = []
        for x in subcon.index:
            tmpnames = [x for x in ind_to_gene.get(x, ["None"]) if "Gm" not in x]
            if len(tmpnames) == 0:
                tmpnames = ['None']
            gene_name = tmpnames[0]           
            genenames.append(gene_name)
        subcon.index = genenames
        subcon.columns = genenames

        osubmat = subcon.values.astype(float)
        ax.matshow(osubmat, cmap=cmap, vmin=-1, vmax=1)
        # ax.set_xticks([])
        # ax.set_yticks([])
        inds_to_label = np.where(subcon.index.isin(['Ikzf2', 'Ctla4', 'Tlr5', 
                                                    'Ccdc93', 'Dusp23', 'Idh1', 
                                                    'Gsta3', 'Bcl2',
                                                    'Arl4c', 'Cxcr4', 'Stat4', 'Fasl', 'Lct',
                                                    'Il18rap',
                                                    ]))[0]
        ax.tick_params(axis='y', labelright=True, labelleft=False, left=False, right=True)
        ax.set_xticks([])
        ax.set_yticks(inds_to_label)
        ax.set_yticklabels(subcon.index[inds_to_label])
        # ax.set_yticklabels(subcon.index)
        ax.set_title(f"Bins with megaloops")
        ax.set_xlabel(f'chr{chrom}: n={len(osubmat)}')
        plt.tight_layout()
        us = [0, 1, -1]
        ax.legend([mpatches.Patch(color=cmap((b+1)/2)) for b in us],
               [labeldict[b] for b in us], bbox_to_anchor=(1, 1), loc='upper left', fontsize=12)
        axs[1].remove()
        fig.savefig(f'./plots/figure2/intra_hubs_chr{chrom}.svg')
        # break
    cluster_df['Inds'] = inds_to_keep
    cluster_df['Clusts'] = clusts_of_inds
    full_cluster_df['Inds'] = full_cluster_df['Inds'].astype(int)
    full_cluster_df['Clust'] = full_cluster_df['Clust'].astype(int)

    return cluster_df, full_cluster_df

import time
def create_cluster_mat_NEW(valdf, clust1s, clust2s, all_intra_connections, int_to_chip_clust):
    n = (np.max(valdf['ChIP_Clust']))+1; 
    jointloop_mat = np.zeros((n, n))
    for n_loops, j, clust1 in zip(valdf['n_MegaLoops'], valdf['Inds'], clust1s):
        if n_loops == 0:
            continue
        indsoi = scipy.sparse.find(all_intra_connections[j, :])[1]
        for i in indsoi:
            clust2 = int(valdf.loc[i]['ChIP_Clust'])
            jointloop_mat[clust1, clust2] += 1
    df = pd.DataFrame(jointloop_mat)
    df.columns = pd.Series(df.columns).apply(lambda x: int_to_chip_clust[x])
    df.index = df.columns    
    return df

def create_cluster_mat(valdf, clust1s, clust2s, all_intra_connections, int_to_chip_clust):
    n = (np.max(valdf['ChIP_Clust']))+1; 
    jointloop_mat = np.zeros((n, n))
    for j, clust1 in zip(valdf['Inds'], clust1s):
        for i, clust2 in zip(valdf['Inds'], clust1s):
            if all_intra_connections[i, j]:
                jointloop_mat[clust1, clust2] += 1
    df = pd.DataFrame(jointloop_mat)
    df.columns = pd.Series(df.columns).apply(lambda x: int_to_chip_clust[x])
    df.index = df.columns    
    return df
