import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns
import scipy
import sklearn
from sklearn.decomposition import PCA
from make_figure3 import *
import hichip_volcanoes
from hichip_volcanoes import *
from tad_functions import *
from scipy.stats import rankdata 


arr = np.asarray
def plot_focal_contacts(inter_and_intra_connections_treg, inter_and_intra_connections_tcon,
                        inter_and_intra_connections_all, all_ind_to_region,
                        goodinds, o, clusts, le, label, filt=False):    
    np.random.seed(2)
    submat = inter_and_intra_connections_treg[goodinds, :][:, goodinds].copy().astype(float)
    submat = submat[o, :][:, o]
    if filt == True:
        submat = scipy.ndimage.gaussian_filter(submat, sigma=1)
    n_clusts = len(np.unique(clusts)); print("N_clusts:", n_clusts)
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    ax = axs[0]
    ax.matshow(submat, cmap=cm.bwr, vmin=-2, vmax=2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"Treg Focal Contacts")
    newax = ax.inset_axes(transform=ax.transAxes,
           bounds = (1.025, 0, .05, 1))
    newax.matshow(clusts[o][:, None], cmap=cm.tab20, vmin=0, vmax=n_clusts, aspect='auto')
    newax.set_title("Hub-like Clustering", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)

    newax.set_xticks([])
    newax.set_yticks([])

    ax = axs[1]
    
    focal_submat = inter_and_intra_connections_tcon[goodinds, :][:, goodinds].copy().astype(float)
    focal_submat = focal_submat[o, :][:, o]
    if filt == True:
        focal_submat = scipy.ndimage.gaussian_filter(focal_submat, sigma=1)
    
    ax.matshow(focal_submat, cmap=cm.bwr, vmin=-2, vmax=2)
    ax.set_xticks([])
    ax.set_yticks([])

    newax = ax.inset_axes(transform=ax.transAxes,
           bounds = (1.05, 0, .05, 1))
    newax.matshow(clusts[o][:, None], cmap=cm.tab20, vmin=0, vmax=n_clusts, aspect='auto')
    newax.set_xticks([])
    newax.set_yticks([])
    newax.set_title("Hub-like Clustering", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)

    chromlist = [all_ind_to_region[x] for x in goodinds]
    chromlabels = le.transform(arr(chromlist)[:, 0])

    newax = ax.inset_axes(transform=ax.transAxes,
           bounds = (1.15, 0, .05, 1))
    newax.matshow(chromlabels[o][:, None], cmap=cm.tab20, aspect='auto')    
    
    
    
    newax.set_xticks([])
    newax.set_yticks([])
    newax.set_title("Chromosome", fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)
    ax.set_title(f"Tcon Focal Contacts")
    
    
    ax = axs[2]
    
    focal_submat = inter_and_intra_connections_all[goodinds, :][:, goodinds].copy().astype(float)
    focal_submat = focal_submat[o, :][:, o]
    if filt == True:
        focal_submat = scipy.ndimage.gaussian_filter(focal_submat, sigma=1)
    
    ax.matshow(focal_submat, cmap=cm.bwr, vmin=-2, vmax=2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"All Megaloops")
    fig.suptitle(f'Clustering on: {label}')
    return fig, ax, focal_submat


def plot_focal_contacts_final(inter_and_intra_connections_all, all_ind_to_region,
                        goodinds, o, clusts, le, compartment, SE_count_dict, bw_val_dict, filt=False):    
    np.random.seed(2)
    
    n_clusts = len(np.unique(clusts)); print("N_clusts:", n_clusts)
    fig, axs = plt.subplots(figsize=(6, 6))
    ax = axs
    
    unordered_focal_submat = inter_and_intra_connections_all[goodinds, :][:, goodinds].copy().astype(float).copy()
    focal_submat = unordered_focal_submat[o, :][:, o]
    if filt == True:
        focal_submat = scipy.ndimage.gaussian_filter(focal_submat, sigma=1)
    
    ax.matshow(focal_submat, cmap=cm.bwr, vmin=-2, vmax=2)
    ax.set_xticks([])
    ax.set_yticks([])

    for c, (name, vals) in enumerate(bw_val_dict.items()):
        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (-.56 - c*.08, 0, .05, 1))
        newax.matshow(np.log2(vals[o][:, None]), cmap='bwr', aspect='auto')
        newax.set_xticks([])
        newax.set_yticks([])
        newax.set_title(name, fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)

    for c, (cond, v) in enumerate(SE_count_dict.items()):
        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (-.48+c*.08, 0, .05, 1))
        newax.matshow(v[o][:, None], cmap='bwr', vmin=-1, vmax=1, aspect='auto')
        newax.set_xticks([])
        newax.set_yticks([])
        newax.set_title(cond, fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)

    chromlist = [all_ind_to_region[x] for x in goodinds]
    chromlabels = le.transform(arr(chromlist)[:, 0])

    valdict = {
        'Compartment' : compartment,
        'Clustering' : clusts,
        'Chromosome' : chromlabels,
    }
    print(clusts.max())
    cluster_custom_cmap = sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20)
    chrom_custom_cmap = sns.color_palette('tab20c', n_colors=20) + sns.color_palette('tab20b', n_colors=20)
    kwargdict = {
        'Compartment' : dict(vmin=-1.5, vmax=1.5, cmap='bwr'),
        'Clustering' : dict(cmap=ListedColormap(cluster_custom_cmap), vmin=0, vmax=39),
        'Chromosome' : dict(cmap=ListedColormap(chrom_custom_cmap), vmin=0, vmax=39),
    }

    for c, (cond, v) in enumerate(valdict.items()):
        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (-.24+c*.08, 0, .05, 1))
        if cond == 'Clustering':
            newax.matshow(v[o][:, None],  aspect='auto', **kwargdict[cond])
        else:
            newax.matshow(v[o][:, None],  aspect='auto', **kwargdict[cond])
        newax.set_xticks([])
        newax.set_yticks([])
        newax.set_title(cond, fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.98)
    
    ax.set_title(f"All Megaloops")
    # fig.suptitle(f'Clustering on: {label}')
    return fig, ax, unordered_focal_submat

from matplotlib.colors import LinearSegmentedColormap
def plot_focal_contacts_submat_final(
                        inter_and_intra_connections_treg, 
                        inter_and_intra_connections_tcon,
                        all_ind_to_region, goodinds, o, og_clusts, clusts, le, 
                        compartment, SE_count_dict, bw_val_dict, clustsoi, filt=False, genes_to_plot=[], ind_to_gene=[],
                        plot_ets1_se=True):
    np.random.seed(2)
    print(bw_val_dict.keys())

    cmap = LinearSegmentedColormap.from_list(
        name='test', 
        colors=['blue', 'white', 'lightgray','white','purple',]
    )    
    cmap.set_bad(color='white')
    n_clusts = min(len(np.unique(clusts)), 20); print("N_clusts:", n_clusts)
    fig, axs = plt.subplots(figsize=(8, 8))
    ax = axs
    
    subinds = np.isin(clusts, clustsoi)[o]
    # print("WE ARE HEREEEE")
    # print('YOOO', np.sum(goodinds[o][subinds] == 5145))
    
    focal_submat = inter_and_intra_connections_treg[goodinds, :][:, goodinds].copy().astype(float)
    focal_submat2 = inter_and_intra_connections_tcon[goodinds, :][:, goodinds].copy().astype(float)
    final_focal_submat = np.zeros_like(focal_submat)

    final_focal_submat[(focal_submat == 0) & (focal_submat2 == 0)] = np.nan
    final_focal_submat[(focal_submat == 1) & (focal_submat2 == 1)] = 0
    final_focal_submat[(focal_submat == 1) & (focal_submat2 == 0)] = 1
    final_focal_submat[(focal_submat == 0) & (focal_submat2 == 1)] = -1


    focal_submat = final_focal_submat[o, :][:, o][subinds, :][:, subinds]
    
    if filt == True:
        focal_submat = scipy.ndimage.gaussian_filter(focal_submat, sigma=1)
    
    sns.heatmap(focal_submat, cmap=cmap, vmin = -1.2, vmax = 1.2, ax=ax, cbar=False, rasterized=True,
                # interpolation='None'
                )
    ax.set_xticks([])
    ax.set_yticks([])

    chromlist = [all_ind_to_region[x] for x in goodinds]
    chromlabels = le.transform(arr(chromlist)[:, 0])

    bw_width = .02
    for c, (name, vals) in enumerate(bw_val_dict.items()):
        name = name.replace("Treg ", "")
        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (1+.56 - c*(bw_width*1.5), 0, bw_width, 1))
        newax.matshow(scipy.stats.zscore(vals)[o][subinds][:, None], cmap='coolwarm', aspect='auto',
                        vmin=-3, vmax=3)

        newax.set_xticks([])
        newax.set_yticks([])
        newax.set_title(name, fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.99)
        for spine in newax.spines.values():
            spine.set_visible(False)

    for c, (cond, v) in enumerate(SE_count_dict.items()):
        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (1+.24+c*(bw_width*1.5), 0, bw_width, 1))
        newax.matshow(v[o][subinds][:, None], cmap='gray_r', vmin=0, vmax=1, aspect='auto')
        newax.set_xticks([])
        newax.set_yticks([])
        newax.set_title(cond, fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.99)
        for spine in newax.spines.values():
            spine.set_visible(False)
        break

    valdict = {
        'OG Clustering' : og_clusts,
        'Clustering' : rankdata(clusts, method='dense')-1,
        'Chromosome' : chromlabels,
        'Compartment' : compartment,
    }


    cluster_palette = sns.color_palette(['lightgreen', 'green', 'orange'])
    og_cluster_palette = sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20)
    # Create a custom colormap from the palette
    cluster_cmap = ListedColormap(cluster_palette)

    kwargdict = {
        'OG Clustering' : dict(cmap=ListedColormap(og_cluster_palette), vmin=0, vmax=39),
        'Clustering' : dict(cmap=cluster_cmap, vmin=0, vmax=2),
        # dict(vmin=0, vmax=n_clusts, cmap='tab20'),
        'Chromosome' : dict(cmap='tab20'),
        'Compartment' : dict(vmin=-1.5, vmax=1.5, cmap='coolwarm'),
    }

    for c, (cond, v) in enumerate(valdict.items()):
        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (1+.1+c*(bw_width*1.5), 0, bw_width, 1))
        newax.matshow(v[o][subinds][:, None],  aspect='auto', **kwargdict[cond])
        newax.set_xticks([])
        newax.set_yticks([])
        newax.set_title(cond, fontsize=10, rotation=40, ha='left', va='bottom', x=-.1, y=.99)
        for spine in newax.spines.values():
            spine.set_visible(False)

    if plot_ets1_se:
        genes_in_order = [ind_to_gene.get(x, []) for x in goodinds[o][subinds]]
        genes_in_order[int(np.where(goodinds[o][subinds] == 5217)[0])] = ['Ets1-SE']
        gene_index = []
        for c, i in enumerate(genes_to_plot):
            FOUND = False
            for c2, j in enumerate(genes_in_order):
                if i in j:
                    gene_index.append(c2)
                    FOUND = True
                    break
            if FOUND == False:
                print("Could NOT find", i)
                # raise Exception
        gene_index.append(np.where(goodinds[o][subinds] == 5217)[0])
        genes_to_plot.append("Ets1-SE")
        

        genes_to_plot_in_order = pd.Series(gene_index, index = genes_to_plot).sort_values()
        genes_to_plot_in_order = genes_to_plot_in_order.index
        n = len(genes_to_plot_in_order)
        for c, i in enumerate(genes_to_plot_in_order):
            for c2, j in enumerate(genes_in_order):
                if i in j:
                    cluster = (rankdata(clusts[o][subinds], method='dense')-1)[c2]
                    color = cluster_cmap(cluster)

                    width = 1/len(focal_submat)/2
                    location = 1 - (c2 / len(focal_submat) + width)
                    ## Draw an arrow from text to location on heatmap:
                    plt.sca(ax)
                    textx, texty = -.05, 1-(c+.5)/n
                    plt.text(textx, texty, i, transform=ax.transAxes, ha='right', va='center', fontsize=8,
                        )
                    plt.annotate('', xy=(textx, texty), xycoords='axes fraction', xytext=(0, location),
                                arrowprops=dict(arrowstyle="-", color='gray', lw=.5))
                    break

    ax.set_title(f"All Megaloops", fontsize=24)
    # fig.suptitle(f'Clustering on: {label}')
    return fig, ax, focal_submat, clusts[o][subinds]


def make_focal_df(connections, goodinds, o, names, indsoi, filt=True):
    focal_submat = connections[goodinds, :][:, goodinds].copy().astype(float)
    focal_submat = focal_submat[o, :][:, o]
    if filt==True:
        focal_submat = scipy.ndimage.gaussian_filter(focal_submat, sigma=1)
    focal_df = pd.DataFrame(focal_submat)
    focal_df.index = names
    focal_df.columns = names
    focal_df = focal_df.iloc[indsoi, :].iloc[:, indsoi]
    return focal_df

def make_odict_and_clustdict(full_cluster_df, all_inter_connections, parsed_chroms, all_ind_to_region,
                                connect_dict, n_clusters):
    # tmpmat = (all_inter_connections)
    # good_inter_inds = np.where(tmpmat.sum(axis=1)>15)[0]
    # good_intra_inds = list(full_cluster_df['Inds'].values)

    tmpmat = (all_inter_connections)
    good_inter_inds = np.where(tmpmat.sum(axis=1)>15)[0]
    intra_indset = set(full_cluster_df[full_cluster_df['goodClust']==1]['Inds'].values)
    inter_indset = set(good_inter_inds)
    goodinds = np.sort(list(intra_indset.union(inter_indset))) 

    le = sklearn.preprocessing.LabelEncoder()
    chroms = arr([all_ind_to_region[x] for x in goodinds])[:, 0]
    le.fit(parsed_chroms)

    chromlist = [all_ind_to_region[x] for x in goodinds]
    chromlabels = le.transform(arr(chromlist)[:, 0])

    clustdict = {}
    odict = {}
    for i, v in connect_dict.items():
        submat = v[goodinds, :][:, goodinds].copy()
        pca = PCA(n_components=20)
        pca = pca.fit(submat)
        components = pca.transform(submat)
        tmp_o, tmp_clusts = make_order2_and_cluster(components, n_clusters=n_clusters)
        tmp_o = tmp_o[::-1]
        tmp_clusts = rename_clusts_by_order(tmp_clusts, tmp_o)
        clustdict[i] = tmp_clusts
        odict[i] = tmp_o
    return clustdict, odict, le, goodinds, intra_indset, inter_indset
    # submat_tcon = inter_and_intra_connections_tcon[goodinds, :][:, goodinds].copy()
    # submat_treg = inter_and_intra_connections_treg[goodinds, :][:, goodinds].copy()
    # jointmat = np.concatenate([submat_tcon, submat_treg], axis=1)
    # pca = PCA(n_components=40)
    # pca = pca.fit(jointmat)
    # components = pca.transform(jointmat)
    # tmp_o, tmp_clusts = make_order_and_cluster(components, n_clusters=10)
    # clustdict['joint'] = tmp_clusts
    # odict['joint'] = tmp_o


def make_odict_and_clustdict_v2(connect_dict, parsed_chroms, all_ind_to_region,
                                indset, 
                                n_clusters, method='ward', metric='euclidean', n_pcs = 20, random_state=1):
    # tmpmat = (all_inter_connections)
    # good_inter_inds = np.where(tmpmat.sum(axis=1)>15)[0]
    # good_intra_inds = list(full_cluster_df['Inds'].values)

    le = sklearn.preprocessing.LabelEncoder()
    le.fit(parsed_chroms)

    clustdict = {}
    odict = {}
    print("HIIO: random state: ", random_state)
    for i, v in connect_dict.items():
        submat = v[indset, :][:, indset].copy()
        pca = PCA(n_components=n_pcs, random_state=random_state)
        pca = pca.fit(submat)
        components = pca.transform(submat)
        tmp_o, tmp_clusts = make_order2_and_cluster(components, n_clusters=n_clusters, method=method, metric=metric)
        # tmp_o, tmp_clusts = make_order2_and_cluster(submat, n_clusters=n_clusters, method=method, metric=metric)

        tmp_o = tmp_o[::-1]
        tmp_clusts = rename_clusts_by_order(tmp_clusts, tmp_o)
        clustdict[i] = tmp_clusts
        odict[i] = tmp_o
    return clustdict, odict, le, indset


def get_bins_interacting_with_megaloop_hub(deseq_effect_mat, deseq_effect_mat_inter, 
                                            deseq_lfc_mat, deseq_lfc_mat_inter,
                                            goodinds, clusts, ):
    intra_submat = deseq_effect_mat.copy()
    intra_submat[np.isnan(intra_submat)] = 0
    inter_submat = deseq_effect_mat_inter.copy()
    inter_submat[np.isnan(inter_submat)] = 0
    effect_mat = (intra_submat+inter_submat)[:, goodinds]
    abs_effects = []
    for i in range(len(effect_mat)):
        tmpeffects = []
        for u in np.unique(clusts):
            tmpeffects.append(np.nanmean(np.abs(effect_mat[i, clusts==u])))
        abs_effects.append(tmpeffects)
        
    eff = arr(abs_effects)
    good_inds2 = np.where(eff.sum(axis=1) > 1)[0]

    l = good_inds2[(eff[good_inds2, :] > .75).any(axis=1)]
    intra_submat = deseq_lfc_mat.copy()
    intra_submat[np.isnan(intra_submat)] = 0
    inter_submat = deseq_lfc_mat_inter.copy()
    inter_submat[np.isnan(inter_submat)] = 0
    full_effect_mat = (intra_submat+inter_submat)

    newclusts = []
    for i in l:
        place = np.where(goodinds==i)[0]
        if len(place) > 0:
            newclusts.append(int(clusts[place]))
        else:
            newclusts.append(np.inf)
    newclusts = arr(newclusts)	
    return full_effect_mat, l, good_inds2, newclusts


def plot_bins_interacting_with_megaloop_hub(full_effect_mat, l, goodinds, 
                                    le, all_ind_to_region, ind_to_gene, inter_and_intra_connections,
                                    clusts, o, newclusts, plot=True):

    subm2 = full_effect_mat[l, :][:, goodinds][:, o]

    np.random.seed(1)
    pca = PCA(n_components=10)
    pca = pca.fit(subm2)
    components = pca.transform(subm2)
    o2, clusts2 = make_order_and_cluster(components, n_clusters=5)
    clusts2 = rename_clusts_by_order(clusts2, o2)

    chromlabels = le.transform(arr([all_ind_to_region[x] for x in goodinds])[:, 0])
    chromlabels2 = le.transform(arr([all_ind_to_region[x] for x in l])[:, 0])

    pca = PCA(n_components=10)
    pca = pca.fit(subm2.T)
    components = pca.transform(subm2.T)
    o2T, _ = make_order_and_cluster(components, n_clusters=10)
    genes_to_label = [   'Lrrc32', 'Lypd6b', 'Pde3b', 'Satb1', 'Itgb8', 
        'Pard3b', 'Hdac9', 'Cxcl3', 'Mdfic', 'Socs2', 'Cd80', 'Tnfrsf8', 'Pdlim4', 'Csf2',
        'Il2ra', 'Dpt', 'Il1rl1', 'Ctla4', 'Il9r']
    lfc_with_cluster_df = pd.DataFrame(subm2)
    lfc_with_cluster_df.index = [get_name(x, ind_to_gene) for x in l]
    lfc_with_cluster_df = lfc_with_cluster_df.iloc[o2, :]

    inds_to_label = np.where(lfc_with_cluster_df.index.isin(genes_to_label))[0]

    if plot == True:
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
        ax = axs[0]
        cax = ax.inset_axes(transform=ax.transAxes,
            bounds = (1, 0, .05, 1))
        sns.heatmap(lfc_with_cluster_df, cmap=cm.bwr, 
                    vmin=-2, vmax=2, ax=ax, rasterized=True, cbar_ax = cax)
        ax.set_title("LFC (Treg ÷ Tcon)")
        ax.set_xticks([]); ax.set_yticks(inds_to_label); 
        ax.set_yticklabels(lfc_with_cluster_df.index[inds_to_label], fontsize=12); 

        offset = -.30
        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (-.15+offset, 0, .05, 1))
        newax.matshow(newclusts[o2][:, None], cmap=cm.tab20, vmin=0, vmax=10, aspect='auto', )
        newax.set_xticks([]); newax.set_yticks([])
        newax.set_title('Hub-like\n Clustering', fontsize=12, rotation=40, ha='left', x=-.25)

        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (-.3+offset, 0, .05, 1))
        newax.matshow(clusts2[o2][:, None], cmap=cm.Accent, vmin=0, vmax=6, aspect='auto', )
        newax.set_xticks([]); newax.set_yticks([])
        newax.set_title('LFC Clustering', fontsize=12, rotation=40, ha='left', x=-.25)

        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (-.45+offset, 0, .05, 1))
        newax.matshow(chromlabels2[o2][:, None], cmap=cm.tab20, aspect='auto', )
        newax.set_xticks([]); newax.set_yticks([])
        newax.set_title('Chromosome', fontsize=12, rotation=40, ha='left', x=-.25)

        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (0, -.15, 1, .05))
        newax.matshow(clusts[o][None, :], cmap=cm.tab20, vmin=0, vmax=10, aspect='auto', )
        newax.set_xticks([]); newax.set_yticks([])
        newax.set_title('Hub-like Clustering', fontsize=12)

        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (0, -.3, 1, .05))
        newax.matshow(chromlabels[o][None, :], cmap=cm.tab20, aspect='auto', )
        newax.set_xticks([]); newax.set_yticks([])
        newax.set_title('Chr', fontsize=12)

        ax = axs[1]
        focal_submat = inter_and_intra_connections[l, :][:, goodinds][:, o][o2, :]
        sns.heatmap(focal_submat, ax=ax, cbar=False, rasterized=True,
                cmap=cm.bwr, vmin=-2, vmax=2); ax.set_xticks([]); ax.set_yticks([]);
        ax.set_title("Megaloops")

        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (0, -.15, 1, .05))
        newax.matshow(clusts[o][None, :], cmap=cm.tab20, vmin=0, vmax=10, aspect='auto',)
        newax.set_xticks([]); newax.set_yticks([])
        newax.set_title('Hub-like Clustering', fontsize=12)

        newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (0, -.3, 1, .05))
        newax.matshow(chromlabels[o][None, :], cmap=cm.tab20, aspect='auto',)
        newax.set_xticks([]); newax.set_yticks([])
        newax.set_title('Chr', fontsize=12)

        return fig, lfc_with_cluster_df, focal_submat, clusts2, o2
    else:
        return lfc_with_cluster_df, focal_submat, clusts2, o2

from plotting_functions import *
from hic_zscore_functions import make_obs_exp_inter 
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def get_interchromosomal_megaloop_hub_pileup(chrom, parsed_chroms, cooldict_50kb, 
                chrom_to_start, chrom_to_end, DIR = './for_snipping/pileups/', index_resolution = 250_000, 
                fetch_resolution = 50_000, add_chr=False):
    
    stride = index_resolution // fetch_resolution                
    all_pileup_dict = {
        x: defaultdict(list) for x in cooldict_50kb.keys()
    }
    for _, cool_key in enumerate(cooldict_50kb.keys()):
        places = []
        for chrom2 in parsed_chroms:
            if chrom2 == chrom:
                continue
            if add_chr:
                treg_data = cooldict_50kb[cool_key].matrix().fetch('chr' + chrom, 'chr' + chrom2)
                treg_data = make_obs_exp_inter(treg_data, nan_to_zero=False, pc = 1e-5)
            else:
                print(chrom, chrom2)
                treg_data = cooldict_50kb[cool_key].matrix().fetch(chrom, chrom2)
                treg_data = make_obs_exp_inter(treg_data, nan_to_zero=False, pc = 1e-5)

            start =  25
            start2 = 25
            
            end = chrom_to_end[chrom] - chrom_to_start[chrom] - 25
            end2 = chrom_to_end[chrom2] - chrom_to_start[chrom2] - 25

            assert abs(treg_data.shape[0]//5 - (end + 25)) < 2
            assert abs(treg_data.shape[1]//5 - (end2 + 25)) < 2

            partner_inds = range(start2, end2)
            base_inds = range(start, end)
            for c, i in enumerate(base_inds):
                treg_mats = []

                for j in partner_inds:
                    sl1 = slice(stride*i-25, stride*i+25)
                    sl2 = slice(stride*j-25, stride*j+25)
                    treg_mat = treg_data[sl1, :][:, sl2]
                    treg_mats.append(treg_mat)

                all_pileup_dict[cool_key][i + chrom_to_start[chrom]] += treg_mats
                if c == 0:
                    places += list(arr(list(partner_inds)) + chrom_to_start[chrom2])
            print("Done with", chrom2)

    import pickle
    from functools import partial

    def save_data(c_key_tuple, all_pileup_dict, places, DIR):
        c, key = c_key_tuple
        for cool_key in all_pileup_dict.keys():
            pickle.dump({key: np.array(all_pileup_dict[cool_key][key])}, open(DIR + f'/ind={key}_{cool_key}_pileup_dict.pkl', 'wb'))
        pickle.dump(np.array(places), open(DIR + f'/ind={key}_inds_to_compare_with.pkl', 'wb'))
        print(f'Done with {key} : {c} / {len(all_pileup_dict[cool_key])}')

    save_data_partial = partial(save_data, all_pileup_dict=all_pileup_dict, places=places, DIR=DIR)
    with ThreadPoolExecutor(max_workers=80) as executor:
        executor.map(save_data_partial, enumerate(all_pileup_dict[cool_key].keys()))

    # return treg_mats



def get_stat_interchromosomal_megaloop_hub_pileup(chrom, parsed_chroms, cooldict_50kb, 
                chrom_to_start, chrom_to_end, stat_inds, DIR = './for_snipping/pileups/', index_resolution = 250_000, 
                fetch_resolution = 50_000, add_chr=False):
    
    stride = index_resolution // fetch_resolution                
    all_pileup_dict ={
        x: defaultdict(list) for x in cooldict_50kb.keys()
    }
    for _, cool_key in enumerate(cooldict_50kb.keys()):
        places = []
        for chrom2 in parsed_chroms:
            if chrom2 == chrom:
                continue
            if add_chr:
                treg_data = cooldict_50kb[cool_key].matrix().fetch('chr' + chrom, 'chr' + chrom2)
                treg_data = make_obs_exp_inter(treg_data, nan_to_zero=False, pc = 1e-5)
            else:
                treg_data = cooldict_50kb[cool_key].matrix().fetch(chrom, chrom2)
                treg_data = make_obs_exp_inter(treg_data, nan_to_zero=False, pc = 1e-5)

            start =  25
            start2 = 25
            
            end = chrom_to_end[chrom] - chrom_to_start[chrom] - 25
            end2 = chrom_to_end[chrom2] - chrom_to_start[chrom2] - 25

            assert abs(treg_data.shape[0]//5 - (end + 25)) < 2
            assert abs(treg_data.shape[1]//5 - (end2 + 25)) < 2

            partner_inds = range(start2, end2)
            places += list(arr(list(partner_inds)) + chrom_to_start[chrom2])

            for c, stat_region in enumerate(stat_inds):
                if stat_region[0] != chrom:
                    continue
                else:
                    i = int(stat_region[1])//250_000
                    if i < 25 or i > end - 25:
                        continue
                treg_mats = []

                for j in partner_inds:
                    sl1 = slice(stride*i-25, stride*i+25)
                    sl2 = slice(stride*j-25, stride*j+25)
                    treg_mat = treg_data[sl1, :][:, sl2]
                    treg_mats.append(treg_mat)

                all_pileup_dict[cool_key][i + chrom_to_start[chrom]] += treg_mats
            print("Done with", chrom2)

    import pickle
    from functools import partial

    def save_data(c_key_tuple, all_pileup_dict, places, DIR):
        c, key = c_key_tuple
        for cool_key in all_pileup_dict.keys():
            pickle.dump({key: np.array(all_pileup_dict[cool_key][key])}, open(DIR + f'/ind={key}_{cool_key}_pileup_dict.pkl', 'wb'))
        pickle.dump(np.array(places), open(DIR + f'/ind={key}_inds_to_compare_with.pkl', 'wb'))
        print(f'Done with {key} : {c} / {len(all_pileup_dict[cool_key])}')

    save_data_partial = partial(save_data, all_pileup_dict=all_pileup_dict, places=places, DIR=DIR)
    with ThreadPoolExecutor(max_workers=20) as executor:
        executor.map(save_data_partial, enumerate(all_pileup_dict[cool_key].keys()))

    # return treg_mats



from copy import deepcopy
from hic_zscore_functions import make_obs_exp_intra, make_obs_exp_inter
def normalize_raw_intra(submat, nan_to_zero=True, pc=1e-8, log=True):
    f, exp = make_obs_exp_intra(submat, nan_to_zero=nan_to_zero, pc=pc, log=log)
    return f

def normalize_raw_inter(submat, nan_to_zero=True, pc=1e-8, log=True):
    f = make_obs_exp_inter(submat, nan_to_zero=nan_to_zero, pc=pc, log=log)
    return f

def make_full_mat(intra_raw, inter_raw):
    no_nan_intra = deepcopy(intra_raw)
    no_nan_inter = deepcopy(inter_raw)
    no_nan_intra[np.isnan(no_nan_intra)] = 0
    no_nan_inter[np.isnan(no_nan_inter)] = 0

    full_zscore_mat = np.zeros(intra_raw.shape)
    full_zscore_mat += no_nan_inter
    full_zscore_mat += no_nan_intra
    
    full_zscore_mat[np.isnan(intra_raw) & np.isnan(inter_raw)] = np.nan
    return full_zscore_mat

def initialize_interchrom_oe_mats(treg_inter_raw, parsed_chroms, inds_to_region):
    sep_oe_mat_treg = np.zeros(treg_inter_raw.shape)
    tot1 = 0
    for i in parsed_chroms:
        tot2 = 0
        n1 = len(inds_to_region[i])
        for j in parsed_chroms:
            n2 = len(inds_to_region[j])
            if j != i:
                m = deepcopy(treg_inter_raw[tot1:tot1+n1, tot2:tot2+n2])
                f = normalize_raw_inter(m, nan_to_zero=False)
                sep_oe_mat_treg[tot1:tot1+n1, tot2:tot2+n2] = f
            tot2+= n2            
        tot1+=n1
    return sep_oe_mat_treg


def get_cluster_specific_megaloop_enrichments_by_counting(treg_megaloop_df, tcon_megaloop_df, inds_with_many_chroms, clustdict, o):
    megaloop_frequency_pval_df = pd.DataFrame()
    megaloop_frequency_enrichment_df = pd.DataFrame()
    for cluster in inds_with_many_chroms:
        indsoi = clustdict['all'][o] == cluster
        if indsoi.sum() == 0:
            continue
        treg_megaloops_in_cluster = treg_megaloop_df.loc[:, indsoi]
        tcon_megaloops_in_cluster = tcon_megaloop_df.loc[:, indsoi]
        
        x = treg_megaloops_in_cluster.sum(axis=1)
        y = tcon_megaloops_in_cluster.sum(axis=1)
        pvals = []
        enrichments = []
        for c in range(len(x)):
            a = x.values[c]+1
            b = y.values[c]+1
            anot = len(x)-a-1
            bnot = len(y)-b-1
            enrichment, pval = scipy.stats.fisher_exact([[a, anot], [b, bnot]])
            pvals.append(pval)
            enrichments.append(enrichment)
        pvals = pd.Series(pvals, index = x.index)
        enrichments = pd.Series(enrichments, index = x.index)
        megaloop_frequency_pval_df[cluster] = pvals
        megaloop_frequency_enrichment_df[cluster] = enrichments
    return megaloop_frequency_enrichment_df, megaloop_frequency_pval_df

def make_rna_cdfs(cluster_hic_enrichment_df, cluster_hic_pval_df, merged_inds_to_subset, ind_to_gene, 
                    peak_gene_dict, columns_to_names, row_colors_dict, pco=.001):
    n = len(merged_inds_to_subset)
    fig, axs = init_subplots_exact(n*len(peak_gene_dict), 2, sharey=True, sharex=True,
                                    fgsz=(40*mm, 40*mm), dpi=300, space=1.4)
    
    for c1, (peak_gene_name, peak_gene_df) in enumerate(peak_gene_dict.items()):
        for c2, cluster in enumerate(merged_inds_to_subset):
            # c = c2*2 + c1
            c = c1*len(merged_inds_to_subset) + c2
            bins_with_diff_signal = (cluster_hic_enrichment_df.apply(np.sign)*(cluster_hic_pval_df < pco)).loc[:, columns_to_names[cluster]]
            bin_dict = {
                'NS' : bins_with_diff_signal.index[bins_with_diff_signal==0],
                'Treg Up' : bins_with_diff_signal.index[bins_with_diff_signal==1],
                'Treg Down' : bins_with_diff_signal.index[bins_with_diff_signal==-1],
            }
            cdict = {
                'NS' : 'black',
                'Treg Up' : 'red',
                'Treg Down' : 'blue',
            }

            geneset = set(peak_gene_df.index)

            base_genes = [] 
            for x in bin_dict['NS']:
                base_genes += ind_to_gene.get(x, [])
            base_genes = list(set(base_genes).intersection(geneset))
            baseline_vals = peak_gene_df.loc[base_genes, 'thickStart']
            all_val_dict = {}
            for name, indsoi in bin_dict.items():
                genes = []
                for x in indsoi:
                    genes.extend(ind_to_gene.get(x, []))
                genes = list(set(genes).intersection(geneset))
                vals = peak_gene_df.loc[genes, 'thickStart']
                v = np.ravel(vals)
                ypoints = np.linspace(0, 1, len(v))
                o = np.argsort(v)
                axs[c].plot(v[o], ypoints, label=f'{name}',
                        color=cdict[name], linewidth=1)
                all_val_dict[name] = v
            add_pval_to_plot(all_val_dict['NS'], all_val_dict['Treg Up'], label='Treg Up:', ax=axs[c], test=scipy.stats.ranksums)
            add_pval_to_plot(all_val_dict['NS'], all_val_dict['Treg Down'], label='Treg Down:', ax=axs[c], yloc=.9, test=scipy.stats.ranksums)
            axs[c].legend(frameon=False, fontsize=10, bbox_to_anchor=(.5, 0), loc='lower left')
            p = '\n'
            axs[c].set_title(f'{columns_to_names[cluster].replace(p, " ")}', color=row_colors_dict[columns_to_names[cluster]], fontsize=18)
            axs[c].grid(False)
            axs[c].set_xlim(-2, 2)
            for spine in axs[c].spines:
                if (spine != 'bottom') and (spine != 'left'):
                    axs[c].spines[spine].set_visible(False)
            axs[c].set_xlabel("RNA LFC")    
            add_xaxis_labels('Tcon', 'Treg', axs[c], fontsize=6)
            if c%n == 0:
                axs[c].set_ylabel(f"{peak_gene_name} CDF",)
    # fig.suptitle("RNA-seq LFC: genes in bins with ∆ in Hi-C w/ cluster", va='bottom')
    # plt.tight_layout()
    return fig


def make_rna_cdf_based_on_different_sets(set_dict, kwarg_dict, peak_gene_dict, ind_to_gene, title='', axs=None):
    n = len(peak_gene_dict)
    
    if axs is None:
        fig, axs = init_subplots_exact(len(peak_gene_dict), 1, sharex=True,
                                        fgsz=(30*mm, 20*mm), dpi=300, space=1.4, as_list=True)
    else:
        fig = None
        assert len(axs) == n
    for c1, (peak_gene_name, peak_gene_df) in enumerate(peak_gene_dict.items()):
            c = c1
        # for name, indset in setdict.items():
            # cdict = {
            #     'NS' : 'black',
            #     'Treg Up' : 'red',
            #     'Treg Down' : 'blue',
            # }

            geneset = set(peak_gene_df.index)
            base_genes = [] 
            for x in set_dict['NS']:
                base_genes += ind_to_gene.get(x, [])
            base_genes = list(set(base_genes).intersection(geneset))
            baseline_vals = peak_gene_df.loc[base_genes, 'thickStart']
            all_val_dict = {}

            pstr_y = 1
            for name, indsoi in set_dict.items():
                genes = []
                for x in indsoi:
                    genes.extend(ind_to_gene.get(x, []))
                genes = list(set(genes).intersection(geneset))
                vals = peak_gene_df.loc[genes, 'thickStart']
                v = np.ravel(vals)
                ypoints = np.linspace(0, 1, len(v))
                o = np.argsort(v)
                kwargs = kwarg_dict[name]
                p = scipy.stats.ranksums(baseline_vals, v)[1]
                pstr = format_pvalue(p)

                if name == 'NS':
                    label = None
                else:
                    label = f'# Genes: {len(v)}'
                axs[c].plot(v[o], ypoints,
                        linewidth=1,
                        label = label,
                        **kwargs
                        )
                        
                if name != 'NS':
                    plt.text(1, .1*pstr_y, pstr, fontsize=6, transform=axs[c].transAxes,
                    ha = 'right',
                    va = 'center', color = kwargs['color'])
                    pstr_y += 1

                all_val_dict[name] = v
            leg = axs[c].legend(frameon=False, loc='upper left', fontsize=5)
            leg._legend_box.align = "left"
            p = '\n'
            if c == 0:
                axs[c].set_title(title)
            axs[c].grid(False)
            axs[c].set_xlim(-2, 2)
            axs[c].set_ylim(0, 1)
            for spine in axs[c].spines:
                if (spine != 'bottom') and (spine != 'left'):
                    axs[c].spines[spine].set_visible(False)
            axs[c].set_xlabel("RNA LFC")    
            add_xaxis_labels('Tcon', 'Treg', axs[c], fontsize=6, y = -.2)
            axs[c].set_ylabel(f"{peak_gene_name} CDF",)
    # fig.suptitle("RNA-seq LFC: genes in bins with ∆ in Hi-C w/ cluster", va='bottom')
    # plt.tight_layout()
    return fig

# def get_cluster_specific_megaloop_enrichments_by_permuting(treg_megaloop_df, tcon_megaloop_df, inds_with_many_chroms, clustdict, o):
#     pvals = []
#     enrichments = []
#     for c in range(len(treg_megaloop_df)):
#         tmpdf = pd.concat([ treg_megaloop_df.loc[c], 
#                             tcon_megaloop_df.loc[c]], axis=1,
#                         )
#         tmpdf = tmpdf[tmpdf.sum(axis=1) > 0]
#         tmpdf.columns=['Treg', 'Tcon']

#         joint_counts = tmpdf.astype(int).value_counts().unstack().reindex(columns=[0, 1], index=[0, 1]).fillna(0)+1
#         enrichment, pval = scipy.stats.fisher_exact(joint_counts)
#         pvals.append(pval)
#         enrichments.append(enrichment)

#     pvals = pd.Series(pvals, index = treg_megaloop_df.index)
#     enrichments = pd.Series(enrichments, index = treg_megaloop_df.index)
#     return enrichments, pvals

import statsmodels
import statsmodels.stats
import statsmodels.stats.multitest
def get_cluster_specific_megaloop_enrichments_by_permuting(treg_inter_interaction_df, tcon_inter_interaction_df, 
                                                           inds_with_many_chroms, clustdict, o):
    cluster_hic_pval_df = pd.DataFrame()
    cluster_hic_enrichment_df = pd.DataFrame()
    for cluster in inds_with_many_chroms:
        indsoi = clustdict['all'][o] == cluster
        if indsoi.sum() == 0:
            continue
        treg_megaloops_in_cluster = treg_inter_interaction_df.loc[:, indsoi]
        tcon_megaloops_in_cluster = tcon_inter_interaction_df.loc[:, indsoi]
        
        pvals = []
        enrichments = []
        for c in range(len(treg_megaloops_in_cluster)):
            a = treg_megaloops_in_cluster.loc[c]
            b = tcon_megaloops_in_cluster.loc[c]

            a = a[~np.isnan(a).values]
            b = b[~np.isnan(b).values]

            goodinds = (a > -2.5) & (b > -2.5)
            a = a[goodinds]
            b = b[goodinds]
            
            _, pval = scipy.stats.ranksums(a, b)
            pvals.append(pval)
            enrichments.append(np.nanmean(a) - np.nanmean(b))
        pvals = np.asarray(pvals)
        idx = ~np.isnan(pvals)
        pvals[idx] = statsmodels.stats.multitest.fdrcorrection(pvals[idx])[1]
        pvals = pd.Series(pvals, index = treg_inter_interaction_df.index)

        enrichments = pd.Series(enrichments, index = treg_inter_interaction_df.index)
        cluster_hic_pval_df[cluster] = pvals
        cluster_hic_enrichment_df[cluster] = enrichments
    cluster_hic_enrichment_df = cluster_hic_enrichment_df.fillna(0)
    return cluster_hic_enrichment_df, cluster_hic_pval_df


def get_lfcs_pvals_from_inds(inds, ind_to_gene, peak_genes_df):
    genes = []
    peak_genes_df = peak_genes_df.set_index('name')
    for i in inds:
        genes += ind_to_gene.get(i, [])
    lfcs = peak_genes_df.loc[genes]['thickStart'].values
    pvals = peak_genes_df.loc[genes]['score'].values
    return lfcs, pvals

mm = 1/2.54/10  # mm in inches

from matplotlib_venn import venn2
def plot_megaloop_degs(deseq_effect_mat, all_intra_megaloops, ind_to_gene, peak_genes, DE_megaloops_cutoff=1, DE_cutoff=4):
    megaloops_up = set(np.where((np.sum((deseq_effect_mat > DE_cutoff) & (all_intra_megaloops), axis=1) > DE_megaloops_cutoff))[0])
    megaloops_down = set(np.where((np.sum((deseq_effect_mat < -DE_cutoff) & (all_intra_megaloops), axis=1) > DE_megaloops_cutoff))[0])
    megaloops_neither = set(np.where(np.sum(all_intra_megaloops, axis=1) > 1)[0]) - megaloops_up - megaloops_down

    genes = []
    for i in megaloops_down - megaloops_up:
        genes += ind_to_gene.get(i, [])
    deg_down = -get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 6).astype(float)

    genes = []
    for i in megaloops_up - megaloops_down:
        genes += ind_to_gene.get(i, [])
    deg_up = -get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 6).astype(float)

    genes = []
    for i in megaloops_up.intersection(megaloops_down):
        genes += ind_to_gene.get(i, [])
    deg_both = -get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 6).astype(float)

    genes = []
    for i in megaloops_neither:
        genes += ind_to_gene.get(i, [])
    deg_neither = -get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 6).astype(float)

    baseline = deg_neither
    conditions = [
        [deg_down, 'red', f'>{DE_megaloops_cutoff} Tcon Megaloops'],
        [deg_up, 'blue', f'>{DE_megaloops_cutoff} Treg Megaloops'],
        [deg_both, 'purple', 'Both'],
        [deg_neither, 'black', 'Neither'],
    ]
    fig1, ax = plt.subplots(figsize = (3, 3))
    for c, (vector, color, label) in enumerate(conditions):
        pval = scipy.stats.ks_2samp(baseline, vector)[1]
        sns.ecdfplot(vector, color=color, ax=ax, label=f'{label}; \nn={len(vector)}; p={format_pvalue(pval)}')
    plt.xlabel("RNA-seq LFC: Treg vs Tcon", fontsize=16)
    plt.title('RNA in megaloop bins', fontsize=18)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Genes in bins with:', fontsize=12)
    plt.xlim([-8, 8])


    fig2, ax = plt.subplots(figsize = (40*mm, 40*mm))
    venn2([megaloops_up, megaloops_down], 
        [f'Bins w/ >{DE_megaloops_cutoff} \nTreg-up megaloops', f'Bins w/ >{DE_megaloops_cutoff} \nTcon-up megaloops ',],
        set_colors=['red', 'blue'], ax=ax)
    plt.title('Differential megaloops')


    genes = []
    for i in megaloops_down - megaloops_up:
        genes += ind_to_gene.get(i, [])
    pval_down = get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 4).astype(float)

    genes = []
    for i in megaloops_up - megaloops_down:
        genes += ind_to_gene.get(i, [])
    pval_up = get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 4).astype(float)

    genes = []
    for i in megaloops_up.intersection(megaloops_down):
        genes += ind_to_gene.get(i, [])
    pval_both = get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 4).astype(float)

    genes = []
    for i in megaloops_neither:
        genes += ind_to_gene.get(i, [])
    pval_neither = get_col(peak_genes.filter(lambda x: x[3] in genes).saveas(), 4).astype(float)

    baseline = deg_neither
    conditions = [
        [pval_down, 'red', f'>{DE_megaloops_cutoff} Tcon Megaloops'],
        [pval_up, 'blue', f'>{DE_megaloops_cutoff} Treg Megaloops'],
        [pval_both, 'purple', 'Both'],
        [pval_neither, 'black', 'Neither'],
    ]
    fig3, ax = plt.subplots(figsize = (3, 3))
    for c, (vector, color, label) in enumerate(conditions):
        # pval = scipy.stats.ks_2samp(baseline, vector)[1]
        ax.scatter(c, np.mean(vector < .05), label=f'{label}; \nn={len(vector)}; p={format_pvalue(pval)}')
    plt.xlabel("RNA-seq LFC: Treg vs Tcon", fontsize=16)
    plt.title('RNA in megaloop bins', fontsize=18)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Genes in bins with:', fontsize=12)
    plt.xlim([-8, 8])

    return fig1, fig2
# def format_pvalue(pval, PCO=.05):
#     if pval < PCO:
#         return f'{pval:.2e}'
#     else:
#         return f'NS'
