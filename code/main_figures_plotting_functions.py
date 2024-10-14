import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
mm = 1 / 2.54 / 10
from plotting_functions import init_subplots_exact

from aux_functions import add_chr_to_bedtool, add_chr_to_anc
import pybedtools as pbt
from collections import defaultdict

def anchor_to_lfc_and_loops(slop_tss):
    all_loops = add_chr_to_bedtool(pbt.BedTool('final_loops/processed_DESEQ/thresh=0/all_loops.csv'))
    all_anchors = add_chr_to_bedtool(pbt.BedTool('final_loops/processed_DESEQ/thresh=0/all_ancs.csv'))
    
    anchor_to_lfcs = defaultdict(list)
    anchor_to_loop = defaultdict(list)
    for row in all_loops:
        L1 = add_chr_to_anc(row[:3])
        L2 = add_chr_to_anc(row[3:6])
        lfc = float(row[7])
        padj = float(row[-1])    
        
        anchor_to_loop[L1].append(list(row))
        anchor_to_loop[L2].append(list(row))
        
        anchor_to_lfcs[L1].append(float(row[7]))
        anchor_to_lfcs[L2].append(float(row[7]))

    gene2anchor = defaultdict(list)
    for i in slop_tss.intersect(all_anchors, wo=True):
        name = i[3]
        anchor = i[4:7]
        gene2anchor[name].append(anchor)

    return dict(anchor_to_loop), dict(anchor_to_lfcs), dict(gene2anchor)


def make_loop_adata_plot(loop_adata, full_pileup_dict):
    indsoi_dict = {
        'Treg' : (loop_adata.obs['loops_Treg_thresh=0'] > 0), 
        'NS' : (loop_adata.obs['loops_ns_thresh=0'] > 0),
        'Tcon' : (loop_adata.obs['loops_Tcon_thresh=0'] > 0), 
    }
    comparisons = {
        'Treg' : 'Rename_Treg_all_no_chrM',
        'Tconv' : 'Rename_Tconv_all_no_chrM',
    }
    fig, axs = init_subplots_exact(6, 3, fgsz=(40*mm/2, 40*mm/2), dpi=200,
                                  sharey=True, space=1.4)
    count = 0
    for _, (key, indsoi) in enumerate(indsoi_dict.items()):
        for c, (name, comp1) in enumerate(comparisons.items()):
            ax = axs[count]
            pileup_1 = full_pileup_dict[comp1][indsoi].copy()
            m1 = np.nanmean(pileup_1, axis=0)
    
            d = 2
            n = len(m1)//2
            v1 = np.nanmean(pileup_1[:, n-d:n+d+1, n-d:n+d+1], axis=(1, 2))
            delta_mat = m1[5:-6, 5:-6]
            n = len(delta_mat)//2
            delta_mat_reshape = delta_mat
            res = 5
            ax.matshow(delta_mat_reshape, vmin=0, vmax=10, cmap='gist_heat_r',
                      extent = [-res*n, res*n, -res*n, res*n], zorder=3)
            
            if count % 2 == 0:
                ax.set_ylabel(key)
            
            if count >= 4:
                ax.set_xticks([-res*n, 0, res*n])
                ax.set_xticklabels([-res*n, 'Anchor', res*n])
            else:
                ax.set_xticks([-res*n, 0, res*n])
                ax.set_xticklabels([])
            ax.set_yticks([-res*n, 0, res*n])
            ax.set_yticklabels([-res*n, 'Anchor', res*n])
            if c ==0:
                ax.tick_params(left=True, labelleft=True)
            else:
                ax.tick_params(labelleft=False)
            ax.tick_params(bottom=True, labelbottom=True, top=False, labeltop=False)
        
            m = len(delta_mat_reshape)
            ax.text(-res*n, -res*(n-1), f'n={len(v1)}', fontsize=6)

            plt.sca(ax)
            plt.xticks(fontsize=6)
            plt.yticks(fontsize=6)
            count += 1
    fig.savefig(f'./plots/paper/fig1/loops.pdf', bbox_inches='tight')


from plotting_functions import *
from volcano_plot import label_points
import itertools

def gene_loop_scatter_plot(gene_dict, anchor_to_lfcs, gene2anchor):
    lfc_df = gene_dict['Resting']

    geneset = set(lfc_df.index)
    rows = []
    for gene, anchors in gene2anchor.items():
        lfcs = np.mean(np.ravel(list(itertools.chain(*[anchor_to_lfcs[tuple(anchor)] for anchor in anchors]))))
        if gene not in geneset:
            continue
        gene_lfc = lfc_df.loc[gene]['thickStart'].mean()
        pval = lfc_df.loc[gene]['score'].mean()
        rows.append([float(gene_lfc), lfcs, gene, pval])

    gene_loop_lfc_df = pd.DataFrame(rows, columns=['gene_lfc', 'loop_lfc', 'gene', 'pval'])
    gene_loop_lfc_df['color'] = (gene_loop_lfc_df['pval'] < .05) * np.sign(gene_loop_lfc_df['gene_lfc'])
    x, y = (gene_loop_lfc_df['gene_lfc'], 
                gene_loop_lfc_df['loop_lfc'])
    fig, ax = init_subplots_exact(1, 1, dpi=200, fgsz=(60*mm, 60*mm))
    for u in [0, 1, -1]:
        idx = gene_loop_lfc_df['color'] == u
        plt.scatter(x[idx], y[idx], s=4, zorder=3, c = gene_loop_lfc_df['color'][idx], cmap = 'coolwarm', vmin=-2, vmax=2)


    plt.axvline(0, color='lightgray', linestyle='--')
    plt.axhline(0, color='lightgray', linestyle='--')
    plt.grid(False)
    plt.xlabel("RNA-seq LFC")
    plt.ylabel("Hi-C Loop LFC")
    plt.title("Differential expression vs Looping")
    add_xaxis_labels('Tcon', 'Treg', plt.gca(), fontsize=6)
    add_yaxis_labels('Tcon', 'Treg', plt.gca(), fontsize=6)
    plt.xticks([-6, -3, 0, 3, 6, 9])
    plt.xlim(-6, 9.5)
    plt.yticks([-2, -1, 0, 1, 2])
    plt.ylim(-2.5, 2.5)

    genes_to_label = ['Lrrc32', 'Ikzf2', 'Foxp3', 'Il2ra', 'Tgfbr3', 'Pde3b', 'Enc1', 'Satb1', 'Dapl1', 'Themis',
                    'Ctla4']
    idx = np.where(gene_loop_lfc_df['gene'].isin(genes_to_label))[0]
    label_points(x, y, idx, gene_loop_lfc_df['gene'][idx].values, fontsize=8)
    add_pearson_to_plot(x, y, plt.gca(), test=scipy.stats.pearsonr, yloc=.98)
    fig.savefig('./plots/paper/fig1/differential_expression_looping.pdf')



import matplotlib.patches as patches
def draw_curve(x0, x1, h, steps=500):
    xs = np.linspace(x0, x1, steps, endpoint=True)
    y_input = ((xs-x0)/(x1-x0)-.5)*2
    ys = np.sqrt(1-np.power(y_input, 2))*(h)
    return xs, ys

def ikzf2_arc_plot(deseq_effect_mat, bw_val_df_all_250kb, my_treg_comp):
    fig, axs = init_subplots_exact(5, 5, fgsz=(140*mm, 10*mm), dpi = 150, yspace=1.4)
    i = np.where(np.abs(deseq_effect_mat[278, :782]) > 4)[0]
    for i1 in i:
        j1 = 278
        i1, j1 = sorted([i1, j1])
        zscore = deseq_effect_mat[i1, j1]
        x, y = draw_curve(i1*250_000, j1*250_000, zscore*5_250_000)
        alpha = np.clip(zscore, -40, 40)
        alpha = np.abs(alpha)/10*.3

        if deseq_effect_mat[i1, j1] > 4:
            ax = axs[0]
            
            ax.plot(x-278*250_000, y+250_000*400, alpha=alpha, color='red', zorder=2)
            if i1 >= 278:
                x0, x1 = x[-1] - 278*250_000-125_000, x[-1]+125_000 - 278*250_000
            else:
                x0, x1 = x[0] - 278*250_000-125_000, x[0]+125_000 - 278*250_000

            for ax in [axs[3]]:
                ax.axvspan(x0, x1, color='lightgray', linewidth=0, ymin=0, ymax=1.5, alpha=.5, zorder=-1)
        elif deseq_effect_mat[i1, j1] < -4:
            ax = axs[1]
            ax.plot(x-278*250_000, -y-250_000*400, alpha=alpha, color='blue', zorder=1)

            if i1 >= 278:
                x0, x1 = x[-1] - 278*250_000-125_000, x[-1]+125_000 - 278*250_000
            else:
                x0, x1 = x[0] - 278*250_000-125_000, x[0]+125_000 - 278*250_000
            for ax in [axs[4]]:
                ax.axvspan(x0, x1, color='lightgray', linewidth=0, ymin=0, ymax=1.5, alpha=.5, zorder=-1)
        # else:
            # ax.axvspan(x[0] - 278*250_000, x[0]+250_000 - 278*250_000, color='black', linewidth=0,)

    axs[0].plot([], [], color='red', zorder=2, label='Treg-Specific')
    axs[1].plot([], [], color='blue', zorder=1, label='Tcon-Specific')
    # ax.plot([0], [0], color='gray', zorder=-1, label='NS')

    ax = axs[2]
    ys_treg = np.ravel(bw_val_df_all_250kb['Treg H3K27ac'][:782])*.5
    ys_treg[ys_treg < .1] = 0
    xs = np.linspace(0, 782, 782)*250_000
    d = 0
    ax.fill_between(xs-278*250_000, [d]*len(ys_treg), ys_treg, label='Treg H3K27ac', alpha=1, linewidth=0, facecolor='salmon')
    ax.set_ylim([0, 2])


    ax = axs[3]
    # ys_treg = np.log2(np.ravel(bw_val_dict_all['CD4 H3K9me3'][:782])+1)
    # ax.fill_between(xs-278*250_000, [-d]*len(ys_treg), -d-ys_treg*250_000*100, label='Treg H3K27ac', alpha=.5, linewidth=0, facecolor='blue')

    ys_treg = np.ravel(my_treg_comp[:782]).copy()
    ys_treg[ys_treg < 0 ] = 0
    ax.fill_between(xs-278*250_000, [0]*len(ys_treg), ys_treg, label='A comp', alpha=1, linewidth=0, facecolor='salmon')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(6, 6), )
    ax.set_yticks([])
    ax.grid(False)
    ax.set_ylim([0, 1.5])

    ax = axs[4]
    ys_treg = np.ravel(my_treg_comp[:782]).copy()
    ys_treg[ys_treg > 0 ] = 0
    ax.fill_between(xs-278*250_000, [0]*len(ys_treg), -ys_treg, label='B comp', alpha=1, linewidth=0, facecolor='lightblue')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(6, 6), )
    ax.set_yticks([])
    ax.grid(False)
    ax.set_ylim([0, 1.5])


    for c, ax in enumerate(axs):
        ax.grid(False)
        ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False)
        ax.tick_params(axis='y', left=False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim([-278*250_000, (782-278)*250_000])
        if c == 4:
            ax.spines['bottom'].set_visible(True)
        else:
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(labelleft=False, bottom=False, labelbottom=False)
    fig.savefig('./plots/paper/fig1/ikzf2_arc_plot.pdf', bbox_inches='tight')


## Figure 2
def metadomain_metadomain_size(all_intra_metadomains, dpi = 100):
    x, y = np.where(np.triu(all_intra_metadomains))
    dist = np.abs((x-y)*250_000)/1e6

    fig = plt.figure(figsize=(48*mm, 22*mm), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])

    sns.histplot(dist, color = 'gray', stat='density', cumulative=False, 
                fill=False, element='step',  ax=ax, linewidth=1)
    plt.axvline(np.nanmedian(dist), color='black', linestyle='--')

    plt.xlim([0, 150])
    plt.xlabel('Metadomain size (Mb)')
    plt.ylabel("Fraction")
    plt.title("Metadomain distance")
    # plt.legend(title='Metadomains', loc='upper right')
    fig.savefig('./plots/paper/fig2/metadomain_distance.pdf', bbox_inches = 'tight',  dpi = 300)    




def count_triplets(mat):
    c = 0
    for i in range(len(mat)):
        indsoi = np.where(mat[i, :])[0]
        for j in indsoi:
            for k in indsoi:
                if (mat[j, k] > 0):
                    c += 1
    return c


def plot_metadomain_network_enrichment(all_intra_tcon_metadomains, CHROMS_TO_USE, chrom_to_start, chrom_to_end, chrom_sizes, cutoffs = [0]):
    n = len(cutoffs)
    fig1, ax1 = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm))
    fig2, ax2 = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm))
    axs = [ax1, ax2]
    for c, cutoff in enumerate(cutoffs):
        texts = []
        for chrom in CHROMS_TO_USE:
            s = chrom_to_start[chrom]
            e = chrom_to_end[chrom]
            z = np.linspace(0, 30_000)
            plt.sca(axs[c])
        
            submat = all_intra_tcon_metadomains[s:e, s:e]
            inds_with_metadomains = submat.sum(axis=0) > cutoff
            submat = submat[inds_with_metadomains, :][:, inds_with_metadomains]
        
            triplet_counts = count_triplets(submat)
            random_triplet_counts = []
            for i in range(500):
                permute_mat = np.random.permutation(submat)
                random_triplet_counts.append(count_triplets(permute_mat))
            
            mean_random_triplet_counts = np.mean(random_triplet_counts)
            lower_95 = np.percentile(random_triplet_counts, 2.5)
            upper_95 = np.percentile(random_triplet_counts, 97.5)
            # Scatter plot
            # plt.scatter(mean_random_triplet_counts, triplet_counts, label='Mean Counts', color='black')
            if chrom == '1':
                t = plt.text(mean_random_triplet_counts, triplet_counts*1.03, f'chr{chrom}', fontsize=6)
            
            # Add 95% confidence interval error bars
            plt.errorbar(mean_random_triplet_counts, triplet_counts,
                        xerr=[[mean_random_triplet_counts - lower_95], [upper_95 - mean_random_triplet_counts]],
                        fmt='o', ecolor='black', elinewidth=1, capsize=1, capthick=1, label='95% Confidence Interval', color='black',
                        markersize=2)
            plt.plot(z, z, color='gray', linestyle='--')
            plt.xlabel("Permuted matrix")
            plt.ylabel("Real matrix")
            plt.title("Megaloop triplets")
            
            x = dict(chrom_sizes)[chrom]
            y = triplet_counts/mean_random_triplet_counts
            plt.sca(axs[c + n])
            plt.scatter(x, y, color = 'black', zorder=3, s = 1)
            plt.ylabel("Enrichment")
            if chrom == '1':
                t = plt.text(x, y*1.03, f'chr{chrom}', fontsize = 6)
            plt.ylim([0, 4.5])
            plt.yticks([0, 1, 2, 3, 4])
            plt.xlabel("Chromosome size")
            plt.title("Triplet enrichment")
            texts.append(t)
    return fig1, fig2

from plot_pvals import add_stat_annotation_boxplot_no_hue
def make_metadomain_plot(all_intra_treg_metadomains, SE_count):
    metadomain_count = all_intra_treg_metadomains.sum(axis=0)
    metadomain_SE_count_df = pd.DataFrame([SE_count, metadomain_count]).T  #.value_counts([0, 1]).unstack()
    metadomain_SE_count_df[0] = metadomain_SE_count_df[0].astype(int)
    metadomain_SE_count_df[metadomain_SE_count_df[0] >= 3] = 2
    
    colors = ['lightgray', 'violet', 'purple']
    fig = plt.figure(figsize=(40*mm, 40*mm), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    sns.boxplot(data = metadomain_SE_count_df, x = 0, y = 1, showfliers=False, ax=ax, zorder=3,
               palette=colors)
    # for x, count in t[0].value_counts().items():
    #     plt.text(x, 30, f'n={count}', fontsize=6, ha='center')
    
    plt.ylim([-1, 35])

    # Create labeled rectangle patches for each color
    labels = [f'n={x}' for x in metadomain_SE_count_df[0].value_counts()]
    patches_list = [patches.Patch(color=color, label=label) for color, label in zip(colors, labels)]
    # Add legend with these patches to the top left corner of the axis
    ax.legend(handles=patches_list, loc='lower left', frameon=False, bbox_to_anchor=(0, .3))
    
    plt.xlabel('# SEs in bin')
    plt.ylabel('Metadomains in bin')
    plt.title("Metadomains by SEs")
    add_stat_annotation_boxplot_no_hue(plt.gca(), metadomain_SE_count_df, 0, 1, [0, 1, 2], [[0, 1], [0, 2]],
                                  ymax = 25, log=True,
                                  h = 1)
    ax.set_axisbelow(True)
    return fig


### Figure 3

import collections
import pycircos
from pycircos import *

from tqdm import tqdm
import collections
import pycircos
from pycircos import *

def interchrom_circos_plot(my_treg_comp, my_tcon_comp, all_inter_metadomains, all_ind_to_region, chrom_to_start, chrom_to_end, chromsizes, bw_val_df_all_250kb,
                           parsed_chroms):
    arcdata_dict = collections.defaultdict(dict)
    for c, region in enumerate(all_ind_to_region):
        chrom, s, e = region
        width = e-s
        chrom = 'chr' + chrom
        if chrom not in arcdata_dict:
            arcdata_dict[chrom]["positions"] = []
            arcdata_dict[chrom]["widths"]    = [] 
            arcdata_dict[chrom]["values"]    = [] 
        arcdata_dict[chrom]["positions"].append(s) 
        arcdata_dict[chrom]["widths"].append(width)
        delta_comp = my_treg_comp[c] - my_tcon_comp[c]
        if np.isnan(delta_comp):
            delta_comp = 0
        rescaled_colorval = (delta_comp+1.5)/3
        rescaled_colorval = max(0, rescaled_colorval)
        rescaled_colorval = min(1, rescaled_colorval)
    #     color = (colmap(rescaled_colorval))
        arcdata_dict[chrom]["values"].append(rescaled_colorval)

    #Set chromosomes
    # import patchworklib as pw 

    Garc    = pycircos.Garc
    Gcircle = pycircos.Gcircle
    circle = Gcircle() 
    for chrom in list(chrom_to_start):
        size = chromsizes[chrom]
        if 'M' in chrom or 'Y' in chrom or 'X' in chrom:
            continue
        arc    = Garc(arc_id='chr' + chrom, size=size, interspace=3, 
                    raxis_range=(1000, 1000), labelposition=60, label_visible=True,
                    )
        circle.add_garc(arc) 
    circle.set_garcs()

    for chrom in parsed_chroms:
        if 'M' in chrom or 'Y' in chrom or 'X' in chrom:
            continue

        s, e = chrom_to_start[chrom], chrom_to_end[chrom]
        vals = np.ravel(my_treg_comp[s:e]).copy()
        vals[np.isnan(vals)] = 0

        vals = scipy.ndimage.gaussian_filter1d(vals, sigma=5)
        circle.fillplot('chr' + chrom, data=vals, positions = np.linspace(0, chromsizes[chrom], len(vals)),
                        raxis_range=(780, 980),
                        facecolor='purple', 
                        # rlim=[.06, 1.6],  
                        #linecolor="royalblue", spine=False
                    )

    for chrom in parsed_chroms:
        if 'M' in chrom or 'Y' in chrom or 'X' in chrom:
            continue

        s, e = chrom_to_start[chrom], chrom_to_end[chrom]

        vals = np.ravel(bw_val_df_all_250kb['Treg H3K27ac'][s:e]).copy()
        vals[np.isnan(vals)] = 0
        vals = scipy.ndimage.gaussian_filter1d(vals, sigma = 8)
        vals[vals < .3] = .3
        circle.fillplot('chr' + chrom, data=vals, positions = np.linspace(0, chromsizes[chrom], len(vals)),
                        raxis_range=(560, 760),
                        facecolor='salmon',
                        # rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],  linecolor="royalblue", spine=False
                    )

    arcdata_dict = collections.defaultdict(dict)
    for i, j in tqdm(list(zip(*np.where(np.triu(all_inter_metadomains))))[::10]):
        chrom1, s1, e1 = all_ind_to_region[i]
        chrom2, s2, e2 = all_ind_to_region[j]
        if chrom1 == 'X' or chrom2 == 'X':
            continue

        source = ('chr' + chrom1, s1, e1, 560)
        destination = ('chr' + chrom2, s2, e2, 560)
        circle.chord_plot(source, destination, facecolor='lightgray',
                        edgecolor='None',
                        linewidth = .1, alpha=.1)
    circle.figure.savefig("plots/paper/fig3/circos.pdf", bbox_inches='tight')



from scipy.stats import zscore
from sklearn.cluster import KMeans
import matplotlib as mpl

def intra_inter_metadomain_clustering(all_intra_treg_metadomains, all_inter_treg_metadomains, 
                                      all_intra_tcon_metadomains, all_inter_tcon_metadomains, bw_val_df_all_250kb,
    k = 3):
    bw_val_df_all_250kb = bw_val_df_all_250kb.copy()
    counts_250kb = pd.DataFrame([all_intra_treg_metadomains.sum(axis=1), all_inter_treg_metadomains.sum(axis=1),
    all_intra_tcon_metadomains.sum(axis=1), all_inter_tcon_metadomains.sum(axis=1)],
                    index = ['Intra Treg', 'Inter Treg', 'Intra Tcon', 'Inter Tcon']
                    ).T
    idx = counts_250kb.index[counts_250kb.sum(axis=1) > 0]

    counts_250kb = counts_250kb.loc[idx, ['Intra Treg', 'Intra Tcon', 'Inter Treg', 'Inter Tcon']]

    color_df = pd.DataFrame()
    # [[sns.color_palette('coolwarm', as_cmap=True)((x+1.5)/3)] for x in my_treg_comp[idx]],
    #                        columns = ['compartments'])
    color_df.index = counts_250kb.index

    half_coolwarm = plt.cm.get_cmap('gray_r', 512)
    newcmp = mpl.colors.ListedColormap(half_coolwarm(np.linspace(0, 0.75, 1000)))

    color_df['H3K27ac'] = [newcmp(x/.55) for x in zscore(bw_val_df_all_250kb['Treg H3K27ac'][counts_250kb.index])]
    color_df['H3K9me3'] = [newcmp(x/.55) for x in zscore(bw_val_df_all_250kb['CD4 H3K9me3'][counts_250kb.index])]
    color_df['H3K27me3'] = [newcmp(x/.55) for x in zscore(bw_val_df_all_250kb['Treg H3K27me3'][counts_250kb.index])]

    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(np.log2(counts_250kb+1))


    # Fitting KMeans with the specified number of clusters
    kmeans = KMeans(n_clusters=k, random_state=0)
    kmeans.fit(df_scaled)
    print("# bins:", len(df_scaled))
    # Reorder clusters from lowest number to highest number:
    # Get the order of the clusters
    order = np.argsort(kmeans.cluster_centers_.sum(axis=1))
    # Get the new labels
    labels = np.zeros_like(kmeans.labels_)
    for i, j in enumerate(order):
        labels[kmeans.labels_ == j] = i
    kmeans.labels_ = labels

    # The labels of the clusters for each point
    labels = kmeans.labels_
    o = np.argsort(labels)
    f = sns.clustermap(np.log2(counts_250kb+1).iloc[o], cmap = 'Purples',  vmin = 1,
                    vmax = 6,
                    figsize=(40*mm*1.15, 60*mm*1.15),
                    colors_ratio=.07,
                    row_colors = color_df,
                    row_cluster = False, 
                    col_cluster=False,
                    yticklabels=False,
                    xticklabels=True,
                    zorder = 3,
                    rasterized=False
                    )
    f.fig.dpi = 300
    f.ax_row_colors.grid(False)
    f.ax_heatmap.grid(False)
    return (f, labels, counts_250kb)


def add_legend(divfac):
    plt.scatter([1], [-.6+1], s=40/divfac, color='black')
    plt.scatter([1], [-.6+1.1], s=80/divfac, color='black')
    plt.scatter([1], [-.6+1.2], s=120/divfac, color='black')

    plt.text(1+.1, -.6+1, s=f'n=40', color='black', va='center', ha='left', fontsize=6)
    plt.text(1+.1, -.6+1.1, s=f'n=80', color='black', va='center', ha='left', fontsize=6)
    plt.text(1+.1, -.6+1.2, s=f'n=120', color='black', va='center', ha='left', fontsize=6)
    
    plt.scatter([1], [-.6+.5+.1], s=40/divfac, color=sns.color_palette("Reds", as_cmap=True)(0.03), edgecolor='black', linewidth=.5)
    plt.scatter([1], [-.6+.5+.2], s=40/divfac, color=sns.color_palette("Reds", as_cmap=True)(.46), edgecolor='black', linewidth=.5)
    plt.scatter([1], [-.6+.5+.3], s=40/divfac, color=sns.color_palette("Reds", as_cmap=True)(.66), edgecolor='black', linewidth=.5)
    
    plt.text(1+.1, -.6+.5+.1, s=f'1', color='black', va='center', ha='left', fontsize=6)
    plt.text(1+.1, -.6+.5+.2, s=f'14', color='black', va='center', ha='left', fontsize=6)
    plt.text(1+.1, -.6+.5+.3, s=f'20', color='black', va='center', ha='left', fontsize=6)
    
def draw_network(G, node_colors, node_sizes, pos, ax, key='', draw_selfLoops = False, divfac=4,
                linewidths=.25):
    plt.sca(ax)
    nx.draw(G, node_color = node_colors, 
            pos=pos,
           width=0,
           node_size = node_sizes,
            arrows=False,
           )

    nx.draw_networkx_nodes(G, node_color = node_colors, 
            pos=pos,
           node_size = arr(node_sizes)*1.1,
            edgecolors='black',
            linewidths=linewidths,
           )
    edge_weights = nx.get_edge_attributes(G, 'weight')
    if draw_selfLoops:
        edgelist = list(nx.selfloop_edges(G))
        edgeweights = [edge_weights[x]/divfac*2 for x in edgelist]
        nx.draw_networkx_edges(G, pos, edgelist=edgelist, 
                               arrowstyle='-|>',
                               width=edgeweights,
                               arrows=True,
                               arrowsize = 4,
                              )
    
    edgelist = set(list(nx.edges(G))) - set(list(nx.selfloop_edges(G)))
    edgeweights = [edge_weights[x]/divfac*8 for x in edgelist]
    nx.draw_networkx_edges(G, pos, edgelist=edgelist, width=edgeweights,
                           arrows=False,
                          );
    plt.title(key)

import networkx as nx
import community as community_louvain
def make_networkx_plots(self, all_inter_megaloops):
    seed = 2017
    np.random.seed(seed)
    
    final_focal_submat = self.focal_submat.copy()
    
    m = self.n_megaloop_by_cluster_df.fillna(0).copy()
    m[m < .05] = 0
    n = len(final_focal_submat)
    G = nx.from_numpy_array(m.values)
    pos = nx.spring_layout(G, weight='weight')
    res = 1
    partition = community_louvain.best_partition(G, weight='weight', resolution = res)

    node_alpha = [30+10*(self.n_different_chroms_within_cluster_df.loc[node]>0).sum() for c, node in enumerate(G.nodes())]
    node_alpha = np.asarray(node_alpha)/np.max(node_alpha)
    
    node_colors = [sns.color_palette('Reds', as_cmap=True)((self.n_different_chroms_within_cluster_df.loc[node]>0).sum()/30) for c, node in enumerate(G.nodes())]
    divfac = 6
    node_colors = [(sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20))[node] for c, node in enumerate(G.nodes())]
    node_sizes = [(self.clustdict['all']==node).sum()/divfac for c, node in enumerate(G.nodes())]
    
    node_colors = [sns.color_palette('tab20')[x%20] for x in partition.values()]

    cluster_custom_cmap = sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20)    
    clusters_to_merge = np.unique(self.clustdict['all'][np.isin(self.merged_clustdict['all'],
            self.merged_inds_to_subset)])

    fig, axs = init_subplots_exact(4, 2, dpi=200, fgsz=(30*mm, 30*mm), xspace=1.5)
    for c, key in enumerate(['Treg H3K27ac', 'Treg H3K27me3']):
        node_colors = [np.nanmean(self.bw_val_dict[key][self.clustdict['all']==node]) for c, node in enumerate(G.nodes())]
        if key == 'Treg H3K27ac':
            node_colors = arr(node_colors)
            node_colors[node_colors > 2] = 2
        edge_weights = nx.get_edge_attributes(G, 'weight')
        scaled_edge_weights = [weight * 6 for weight in edge_weights.values()]  # scale factor can be adjusted
        draw_network(G, node_colors, arr(node_sizes)*2, pos, axs[c+2], key=key, linewidths=.75)
    
    
    node_colors = [((self.n_different_chroms_within_cluster_df.loc[node]>0).sum()/30) for c, node in enumerate(G.nodes())]
    node_colors = []
    for x in partition.keys():
        color = cluster_custom_cmap[x%20]
        node_colors.append(color)
    draw_network(G, node_colors, arr(node_sizes)*2, pos, axs[0], key='Clusters', linewidths=.5)
        
        
    node_colors = [sns.color_palette('Reds', as_cmap=True)((self.n_different_chroms_within_cluster_df.loc[node]>0).sum()/30) for c, node in enumerate(G.nodes())]
    draw_network(G, node_colors, arr(node_sizes)*2, pos, axs[1], key='# Chroms', linewidths=.5)
    fig.savefig('./plots/paper/fig3/networkx_plot1.pdf', bbox_inches='tight')


    cluster_custom_cmap = sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20)    
    clusters_to_merge = np.unique(self.clustdict['all'][np.isin(self.merged_clustdict['all'],
            self.merged_inds_to_subset)])
    fig, axs = init_subplots_exact(1, 1, dpi=200, fgsz=(30*mm, 30*mm), xspace=1.5, as_list=True)
    node_colors = []
    for x in partition.keys():
        if x in clusters_to_merge:
            color = cluster_custom_cmap[x%20]
        else:
            color = sns.color_palette('coolwarm', n_colors=5)[2]
        node_colors.append(color)
    draw_network(G, node_colors, arr(node_sizes)*2, pos, axs[0], key='Clusters', linewidths=.5)


    
    j = pd.DataFrame(all_inter_megaloops[self.goodinds, :][:, self.goodinds])
    j['cluster'] = self.clustdict['all']
    j = j.fillna(0).groupby('cluster').sum().T
    j['cluster'] = self.clustdict['all']
    j = j.fillna(0).groupby('cluster').sum().T
    
    
    chromequal = np.equal.outer([self.all_ind_to_region[x][0] for x in self.goodinds], [self.all_ind_to_region[x][0] for x in self.goodinds])
    j2 = pd.DataFrame(~chromequal)
    j2['cluster'] = self.clustdict['all']
    j2 = j2.fillna(0).groupby('cluster').sum().T
    j2['cluster'] = self.clustdict['all']
    j2 = j2.fillna(0).groupby('cluster').sum().T
    
    chromequal = np.equal.outer([self.all_ind_to_region[x][0] for x in self.goodinds], [self.all_ind_to_region[x][0] for x in self.goodinds])
    j3 = pd.DataFrame(~chromequal)
    j3['cluster'] = self.clustdict['all']
    j3 = j3.fillna(0).groupby('cluster').mean().T
    j3['cluster'] = self.clustdict['all']
    j3 = j3.fillna(0).groupby('cluster').mean().T
    
    x = np.diag(j)/(np.diag(j2) + np.diag(j))
    x[np.isnan(x)] = 0
    # entropy = (self.n_different_chroms_within_cluster_df.fillna(0)>0).sum(axis=1)
    entropy = scipy.stats.entropy(self.n_different_chroms_within_cluster_df.fillna(0), axis=1)
    y = entropy

    goodones = (x > .05) & (y > 1)
    colors = sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20)
    fig, axs = init_subplots_exact(1, 1, fgsz=(25*mm, 25*mm), dpi=200, xspace=1.3)
    plt.scatter(x, y , zorder=3, c = colors[:len(x)],  s = 40)
    plt.scatter(x[goodones], y[goodones], zorder=3, c = 'None', edgecolor='black', s = 40,
               label='Hubs')
    plt.axvline(.05, color='black', linestyle='--')
    plt.axhline(2, color='black', linestyle='--')
    patch = matplotlib.patches.Rectangle((.05, 2), 1, 1, color='lightgray', alpha=.6)
    plt.gca().add_patch(patch)

    plt.xlabel("Interchrom. pairs \nw/ megaloop")
    plt.ylabel('Chrom. Entropy')
    # add_yaxis_labels('Only one chrom', 'More chroms', axs, fontsize=4)
    plt.title("Interchrom. hubs")
    plt.legend(bbox_to_anchor=(1, 1), frameon=False, loc = 'upper left')
    # plt.yscale('log')
    fig.savefig('./plots/paper/fig3/networkx_plot2.pdf', bbox_inches='tight')


def chip_boxplot(_250kb_hub_annotations, bw_val_df_all_250kb, my_treg_comp, row_colors, row_colors_dict,
                 hue_order = ['matched_A comp.', 'Constitutive', 'Dynamic', 'Repressive'],
                 palette = ['skyblue'], chip_keys =  ['H3K4me3', 'H3K4me1', 'H3K27ac', 'Smc1a', 'CTCF', 'H3K27me3'],
                 ):
    palette = palette + row_colors
    bw_val_df_all_250kb = bw_val_df_all_250kb.copy()
    treg_chip_keys = ['Treg ' + x for x in chip_keys]
    bw_val_df_all_250kb = bw_val_df_all_250kb[treg_chip_keys].apply(zscore, axis=0)

    v_df = pd.concat([_250kb_hub_annotations, bw_val_df_all_250kb], axis=1)
    v_df['Compartment'] = my_treg_comp

    v_df.loc[(v_df['Hub'] == -1), 'Hub'] = 'Others'
    data = v_df.dropna().melt(id_vars='Hub')
    data['Hub'] = data['Hub']
    data['variable'] = data['variable'].str.replace("Treg ", "")
    data = data[data['Hub']!=-1]

    fig, axs = init_subplots_exact(1, 1, dpi = 200, fgsz=(100*mm, 30*mm), yspace = 1.5, xspace=1.2, sharey=True)
    sns.boxplot(data=data,
                x='variable', y='value', hue='Hub',
                order =  chip_keys,
                hue_order = ['Others'] + hue_order,
                palette = ['lightgray'] + palette,
                # orient='h', 
                showfliers=False, zorder=3,)
    plt.legend(bbox_to_anchor=(0, -.3), loc = 'center left', frameon=False, fontsize=6, ncol=13)
    plt.xlabel("")
    plt.ylabel("ChIP Enrichment")
    plt.title('ChIP Enrichment at Clusters')
    axs.set_axisbelow(True)


    from scipy.stats import ranksums
    from statsmodels.stats.multitest import multipletests
    from matplotlib.lines import Line2D

    sizemap = {
        1e-50 : 12,
        1e-5 : 6,
        .05 : 2,
    }

    def compute_pvalues_and_annotate(ax, data, x, y, hue, baseline_hue, order, hue_order):
        # Dictionary to store p-values and group coordinates
        p_values = {}
        # Iterate over x values and compute p-values
        for i, x_val in enumerate(order):
            baseline_data = data[(data[x] == x_val) & (data[hue] == baseline_hue)][y]

            for j, hue_val in enumerate(hue_order):
                if hue_val != baseline_hue:
                    comparison_data = data[(data[x] == x_val) & (data[hue] == hue_val)][y]
                    _, p_value = ranksums(baseline_data, comparison_data)
                    p_values[(i, j)] = p_value
        # Apply FDR correction
        p_vals = list(p_values.values())
        _, corrected_p_vals, _, _ = multipletests(p_vals, method='fdr_bh')

        # Annotate plot with asterisks for significant comparisons
        for (x_idx, hue_idx), corrected_p_val in zip(p_values.keys(), corrected_p_vals):
            if corrected_p_val < 0.05:  # Check for significance
                # Scale asterisk size by -log10 of the p-value
                for key, val in sizemap.items():
                    if corrected_p_val < key:
                        asterisk_size = val
                        break
                    
                # asterisk_size = max(2, np.sqrt(-np.log10(corrected_p_val) * 10)/3)
                x_coord = x_idx + hue_idx * .18 - 0.26  # Adjust based on number of hue categories

                # Adjust vertical position based on asterisk size
                y_offset = 9# - (asterisk_size / 30)  # Adjust this factor as needed
                ax.scatter(x_coord, -4 + 0*y_offset, s=asterisk_size, color = row_colors_dict[hue_order[hue_idx]],
                        clip_on=False,
                        )
                # Annotate with colored asterisk

        # Custom legend for asterisk sizes
        for c, (key, asterisk_size) in enumerate(sizemap.items()):
            ax.scatter(4, -6-c, s = asterisk_size, clip_on=False, color='black')
            ax.text(4.2, -6-c, f'p<{key}', fontsize = 6, va='center')

        # Adjust ylim to make space for asterisks
        ax.set_ylim(-2.5, top=9)  # Adjust as necessary


    # Compute p-values and annotate

    compute_pvalues_and_annotate(axs, data, 'variable', 'value', 'Hub', 'matched_A comp.',
                                ['H3K4me3', 'H3K4me1', 'H3K27ac', 'Smc1a', 'CTCF', 'H3K27me3'],
                                hue_order)

    # Finalize plot
    plt.legend(bbox_to_anchor=(0, -.3), loc='center left', frameon=False, fontsize=6, ncol=2)
    plt.xlabel("")
    plt.ylabel("ChIP Enrichment")
    plt.title('ChIP Enrichment at Clusters')
    axs.set_axisbelow(True)
    return fig




def compartment_boxplot(_250kb_hub_annotations, bw_val_df_all_250kb, my_treg_comp, row_colors, row_colors_dict,
                 hue_order = ['matched_A comp.', 'Constitutive', 'Dynamic', 'Repressive'],
                 palette = ['skyblue'], chip_keys =  ['H3K4me3', 'H3K4me1', 'H3K27ac', 'Smc1a', 'CTCF', 'H3K27me3'],
                 ):
    palette = palette + row_colors
    bw_val_df_all_250kb = bw_val_df_all_250kb.copy()
    treg_chip_keys = ['Treg ' + x for x in chip_keys]
    bw_val_df_all_250kb = bw_val_df_all_250kb[treg_chip_keys].apply(zscore, axis=0)

    v_df = pd.concat([_250kb_hub_annotations, bw_val_df_all_250kb], axis=1)
    v_df['Compartment'] = my_treg_comp

    v_df.loc[(v_df['Hub'] == -1), 'Hub'] = 'Others'
    data = v_df.dropna().melt(id_vars='Hub')
    data['Hub'] = data['Hub']
    data['variable'] = data['variable'].str.replace("Treg ", "")
    data = data[data['Hub']!=-1]

    fig, axs = init_subplots_exact(1, 1, dpi = 200, fgsz=(20*mm, 30*mm), yspace = 1.5, xspace=1.2, sharey=True)
    sns.boxplot(data=data,
                x='variable', y='value', hue='Hub',
                order =  ['Compartment'],
                hue_order = ['Others'] + hue_order,
                palette = ['lightgray'] + palette,
                # orient='h', 
                showfliers=False, zorder=3,)
    plt.legend(bbox_to_anchor=(0, -.3), loc = 'center left', frameon=False, fontsize=6, ncol=13)
    plt.xlim([-.5, .5])
    axs.set_axisbelow(True)

    from scipy.stats import ranksums
    from statsmodels.stats.multitest import multipletests
    from matplotlib.lines import Line2D

    sizemap = {
        1e-50 : 12,
        1e-5 : 6,
        .05 : 2,
    }

    def compute_pvalues_and_annotate(ax, data, x, y, hue, baseline_hue, order, hue_order):
        # Dictionary to store p-values and group coordinates
        p_values = {}
        # Iterate over x values and compute p-values
        for i, x_val in enumerate(order):
            baseline_data = data[(data[x] == x_val) & (data[hue] == baseline_hue)][y]

            for j, hue_val in enumerate(hue_order):
                if hue_val != baseline_hue:
                    comparison_data = data[(data[x] == x_val) & (data[hue] == hue_val)][y]
                    _, p_value = ranksums(baseline_data, comparison_data)
                    p_values[(i, j)] = p_value
        # Apply FDR correction
        p_vals = list(p_values.values())
        _, corrected_p_vals, _, _ = multipletests(p_vals, method='fdr_bh')

        # Annotate plot with asterisks for significant comparisons
        for (x_idx, hue_idx), corrected_p_val in zip(p_values.keys(), corrected_p_vals):
            if corrected_p_val < 0.05:  # Check for significance
                # Scale asterisk size by -log10 of the p-value
                for key, val in sizemap.items():
                    if corrected_p_val < key:
                        asterisk_size = val
                        break
                    
                # asterisk_size = max(2, np.sqrt(-np.log10(corrected_p_val) * 10)/3)
                x_coord = x_idx + hue_idx * .18 - 0.26  # Adjust based on number of hue categories

                # Adjust vertical position based on asterisk size
                y_offset = 9# - (asterisk_size / 30)  # Adjust this factor as needed
                ax.scatter(x_coord, -4 + 0*y_offset, s=asterisk_size, color = row_colors_dict[hue_order[hue_idx]],
                        clip_on=False,
                        )
                # Annotate with colored asterisk

        # Custom legend for asterisk sizes
        for c, (key, asterisk_size) in enumerate(sizemap.items()):
            ax.scatter(4, -6-c, s = asterisk_size, clip_on=False, color='black')
            ax.text(4.2, -6-c, f'p<{key}', fontsize = 6, va='center')

        # Adjust ylim to make space for asterisks
        ax.set_ylim(-2.5, top=9)  # Adjust as necessary


    # Compute p-values and annotate

    compute_pvalues_and_annotate(axs, data, 'variable', 'value', 'Hub', 'matched_A comp.',
                                ['Compartment'],
                                hue_order)

    # Finalize plot
    plt.legend(bbox_to_anchor=(0, -.3), loc='center left', frameon=False, fontsize=6, ncol=2)
    plt.xlabel("")
    plt.ylabel("Compartment Score")
    plt.title('Compartment Score at Clusters')

    axs.set_axisbelow(True)
    print(1)
    plt.ylim([-2, 2])

    return fig













def rpkm_per_hub(hub_annotations, rpkm_df, row_colors_dict, gene_to_ind, ind_to_gene):
    subgenes = rpkm_df[rpkm_df.index.isin(gene_to_ind.keys())]
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)

    baseline = subgenes

    cluster_valdict = {}
    for cluster in hub_annotations['Hub'].unique():
        if cluster == 'Others':
            continue
        genes = np.unique(list(itertools.chain(*[ind_to_gene.get(x, []) for x in hub_annotations.index[hub_annotations['Hub']==cluster]])))
        genes = pd.Series(genes)
        genes = genes[genes.isin(subgenes.index)]
        rpkm = subgenes.loc[genes]
        p = format_pvalue(scipy.stats.ranksums(rpkm, baseline).pvalue)
        cluster_valdict[cluster] = rpkm
        sns.ecdfplot(rpkm, color = row_colors_dict.get(cluster, 'lightgray'),
                    linewidth = 1,
                    label = f'{cluster}, p={p}')
    
    pval = scipy.stats.ranksums(cluster_valdict['Constitutive'], cluster_valdict['Dynamic'])[1]
    print("Constitutive - Dynamic", pval)    

    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', frameon=False)
    plt.xlabel("Log RPKM (kb)")
    plt.title("Gene expr. in hub")
    plt.xlim([-.2, 4])
    return cluster_valdict




def rpkm_per_hub2(row_colors_dict, col='rTreg'):
    readcounts = {
        'rTreg' : [
            'SRP272473.SRR12264701.resting_Treg.rep1.mapq30.counts.txt', 
            'SRP272473.SRR12264702.resting_Treg.rep2.mapq30.counts.txt', 
            'SRP272473.SRR12264703.resting_Treg.rep3.mapq30.counts.txt',
    ],
        'aTreg' : [
            'SRP272473.SRR12264704.activated_Treg.rep1.mapq30.counts.txt', 
            'SRP272473.SRR12264705.activated_Treg.rep2.mapq30.counts.txt', 
            'SRP272473.SRR12264706.activated_Treg.rep3.mapq30.counts.txt',
    ],
        'aTcon' : [
            'SRP272473.SRR12264698.activated_Tcon.rep1.mapq30.counts.txt',
            'SRP272473.SRR12264699.activated_Tcon.rep3.mapq30.counts.txt',
            'SRP272473.SRR12264700.activated_Tcon.rep4.mapq30.counts.txt',
        ],

        'rTcon' : [
            'SRP272473.SRR12264695.resting_Tcon.rep1.mapq30.counts.txt',
            'SRP272473.SRR12264696.resting_Tcon.rep2.mapq30.counts.txt',
            'SRP272473.SRR12264697.resting_Tcon.rep3.mapq30.counts.txt',
        ]
    }
    rna_pref = './rudensky_scrna/prelim-analysis/bulk_rna_data/counts/'
    rpkm_df = pd.DataFrame()
    for name, files in readcounts.items():
        readcount_data = []
        for file in files:
            d = pd.read_csv(rna_pref + file, sep='\t', skiprows=1)
            d['rpkm'] = d.iloc[:, -1] / d['Length']
            readcount_data.append(d.set_index('Geneid')['rpkm'])
        rpkm_df[name] = pd.concat(readcount_data, axis=1).sum(axis=1)

    all_tss_df = pd.read_csv('./annotations/full_tss_df.csv', index_col=0)
    all_tss_df['hub'] = all_tss_df['hub'].str.replace("compartment", "comp.")

    subgenes = rpkm_df[rpkm_df.index.isin(all_tss_df['gene_name'])][col]
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)

    baseline = subgenes
    cluster_valdict = {}
    for cluster in all_tss_df['hub'].unique():
        genes = all_tss_df['gene_name'][all_tss_df['hub']==cluster]
        genes = pd.Series(genes)
        genes = genes[genes.isin(subgenes.index)]
        rpkm = subgenes.loc[genes]
        p = format_pvalue(scipy.stats.ranksums(rpkm, baseline).pvalue)
        cluster_valdict[cluster] = rpkm
        sns.ecdfplot(rpkm, color = row_colors_dict.get(cluster, 'lightgray'),
                    linewidth = 1,
                    label = f'{cluster}, p={p}')

    pval = scipy.stats.ranksums(cluster_valdict['Constitutive'], cluster_valdict['Dynamic'])[1]
    print("Constitutive - Dynamic", pval)    

    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', frameon=False)
    plt.xlabel("Log RPKM")
    plt.title("Gene expr. in hub")
    plt.xlim([-.2, 4])
    return fig

def make_immgen_datasets(geneLengths, my_tss_df, cluster_key='hub'):
    immgen_data = pd.read_csv('./rudensky_scrna/prelim-analysis/immgen_data/GSE109125_Normalized_Gene_count_table.csv', index_col=0)
    
    geneLengths = geneLengths[~geneLengths.index.duplicated()]
    
    good_genes = immgen_data.index.intersection(geneLengths.index).intersection(my_tss_df['gene_name'])
    immgen_data = immgen_data.loc[good_genes]
    geneLengths = geneLengths.loc[good_genes]
    
    immgen_data = (immgen_data/(immgen_data.sum(axis=0)/1e6)).div(geneLengths, axis=0)*1e3
    
    immgen_metadata = pd.read_csv('./rudensky_scrna/prelim-analysis/immgen_data/immgen_metadata.csv')
    immgen_metadata = immgen_metadata.set_index("SampleName")
    immgen_data = immgen_data.loc[:, immgen_data.columns.isin(immgen_metadata.index)]

    good_genes = immgen_data.index.intersection(geneLengths.index)
    immgen_data = immgen_data.loc[good_genes]
    geneLengths = geneLengths.loc[good_genes]    

    immgen_data = np.log2(immgen_data).apply(zscore, axis=0)

    immgen_data['cluster'] = my_tss_df.set_index('gene_name').loc[immgen_data.index, cluster_key]
    return immgen_data, geneLengths, immgen_metadata

def aggregate_immgen_by_cluster(immgen_data, immgen_metadata, group_lineages=True, cluster_key='cluster'):
    datas = []
    clusters = immgen_data[cluster_key].unique()
    for u in clusters:
        genes = immgen_data.index[immgen_data[cluster_key] == u]
        l = immgen_data.loc[genes].drop(cluster_key, axis=1)
        l = l.mean()
        datas.append(l)
    data = pd.concat(datas, axis=1)
    data.columns = clusters
    data['lineage'] = immgen_metadata.loc[data.index, 'Lineage'].values
    data['celltype'] = immgen_metadata.loc[data.index, 'CellType'].values
    if group_lineages:
        data = data.groupby(['lineage', 'celltype'])[clusters].mean().reset_index()
    return data

import scipy.cluster.hierarchy as sch
def plot_aggregate_immgen_data(immgen_aggregate_data, immgen_gene_expression, row_colors_dict, cluster_key='hub'):
    immgen_lineage_color_dict = dict(zip(
            ['B',        'DC',   'ILC', 'Stem&Prog', 'T.act', 'abT', 'gdT', 'monocyte', 'myeloid', 'stroma'],
            ['purple',   'salmon', 'aqua', 'brown',    'blue', 'blue', 'blue', 'red',     'orange', 'pink']
    ))
    row_colors = immgen_aggregate_data[['lineage']].copy()
    row_colors['lineage'] = row_colors['lineage'].apply(immgen_lineage_color_dict.get)
    row_colors.index = immgen_aggregate_data['celltype']
    clusters = np.unique(immgen_gene_expression[cluster_key])
    values = immgen_aggregate_data[clusters]
    values.index = immgen_aggregate_data['celltype']
    o, _, linkage = make_order_and_cluster_optimal(values, method='ward', metric='euclidean',
                                           n_clusters=10,
                                           return_linkage=True)
    values = values[['Constitutive', 'Dynamic', 'Repressive', 'Others', 'matched_A comp.']]
    values = values.iloc[o]
    g = sns.clustermap(values, zorder=3, cmap='coolwarm', 
                       row_colors=row_colors, 
                       figsize=(4, 4),
                       col_colors=[row_colors_dict[x] for x in values.columns],
                       row_cluster=False,
                       col_cluster=False,
                       metric='euclidean'
                       )
    
    # Manually add the dendrogram
    ax = g.ax_row_dendrogram
    dendro = sch.dendrogram(linkage, ax=ax, orientation='left', no_labels=True, 
                            color_threshold=np.inf,  # Set to inf to prevent coloring
                            link_color_func=lambda k: 'black')
    ax.invert_yaxis()  # Match the order with the heatmap
    g.fig.savefig('./plots/paper/fig3/immgen_hub_score_heatmap.pdf', bbox_inches='tight')    


# Figure 4

# Figure 5
    
def make_differential_chip_plot(pval_df, stat_df, bw_val_df_all_250kb):
    comps = [
    ['Treg H3K27me3', 'Tcon H3K27me3'],
    ['Treg H3K27ac', 'Tcon H3K27ac'],
    ['Treg H3K4me3', 'Tcon H3K4me3'],
    ['Treg H3K4me1', 'Tcon H3K4me1'],
    ]
    delta_chip_df = pd.DataFrame()
    for comp1, comp2 in comps:
        factor = comp1.split(" ")[-1]
        delta_chip_df[factor] = np.ravel(np.log2(1+bw_val_df_all_250kb[comp1]) - np.log2(1+bw_val_df_all_250kb[comp2]))

    pco = .01
    stat_co = 1.5
    for c, chip_factor in enumerate(delta_chip_df.columns[delta_chip_df.columns.str.contains("H3K")]):
        n = len(pval_df.columns)
        fig, axs = init_subplots_exact(n, 1, fgsz=(25*mm, 15*mm), dpi=200, sharey=True, xspace=1.4)
        for _, col in enumerate(pval_df.columns):
            plt.sca(axs[_])
            treg_up = pval_df.index[(pval_df[col] < pco) & (stat_df[col] > stat_co)]
            tcon_up = pval_df.index[(pval_df[col] < pco) & (stat_df[col] < -stat_co)]
            ns_bins = pval_df.index[(pval_df < pco).any(axis=1)]
                
            v1 = zscore(delta_chip_df[chip_factor], nan_policy='omit').loc[treg_up]
            v2 = zscore(delta_chip_df[chip_factor], nan_policy='omit').loc[tcon_up]
            v3 = zscore(delta_chip_df[chip_factor], nan_policy='omit').loc[ns_bins]
            
            name = col
            sns.ecdfplot(v1, color='red', label=f'# Bins: {len(v1)}', linewidth=1,)
            sns.ecdfplot(v2, color='blue', label=f'# Bins: {len(v2)}', linewidth=1,)
            sns.ecdfplot(v3, color='lightgray', linewidth=1,)
    
    
            p1 = scipy.stats.ranksums(v1, v3)[1]
            p2 = scipy.stats.ranksums(v2, v3)[1]
    
            plt.text(1, .1, format_pvalue(p1), fontsize=6, transform=axs[_].transAxes,
                     ha = 'right', va = 'center', color = 'red')
            plt.text(1, .2, format_pvalue(p2), fontsize=6, transform=axs[_].transAxes,
                     ha = 'right', va = 'center', color = 'blue')
    
    
            add_xaxis_labels('Tcon', 'Treg', plt.gca(), fontsize=6)
            plt.legend(frameon=False, fontsize=5, bbox_to_anchor=(-.02, 1.02), loc='upper left')
            
            if _ == 1:
                plt.xlabel(f"Local change in {chip_factor}")
            else:
                plt.xlabel("")
            plt.xlim([-3, 3])
        for ax in axs[1:]:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)
        fig.savefig(f'./plots/paper/fig5/change_in_{chip_factor}_at_diff_megaloop_bins.pdf', bbox_inches='tight')


from make_figure4 import make_rna_cdf_based_on_different_sets
def make_differential_rna_plot(pval_df, stat_df, gene_dict, ind_to_gene, columns_to_names):
    pco = .01
    stat_co = 1.5
    treg_sets = []
    tcon_sets = []
    kwarg_dict = {}
    fig, axs = init_subplots_exact(3, 1, dpi=200, fgsz=(25*mm, 15*mm), xspace=1.4, sharey=True)
    for c, col in enumerate(pval_df):
        treg_up = pval_df.index[(pval_df[col] < pco) & (stat_df[col] > stat_co)]
        tcon_up = pval_df.index[(pval_df[col] < pco) & (stat_df[col] < -stat_co)]
        ns_bins = pval_df.index[((pval_df < pco)).any(axis=1)]
        treg_sets.append(set(treg_up))
        tcon_sets.append(set(tcon_up))

        set_dict = {}
        set_dict[f'Treg Gains'] = treg_up
        set_dict[f'Tcon Gains'] = tcon_up
        set_dict['NS'] = set(ns_bins)
    
        kwarg_dict[f'Treg Gains'] = {'color' : 'red'}
        kwarg_dict[f'Tcon Gains'] = {'color' : 'blue'}
        kwarg_dict[f'NS'] = {'color' : 'lightgray'}
    
        make_rna_cdf_based_on_different_sets(set_dict, kwarg_dict, subset_dict(gene_dict , ['Resting']), ind_to_gene, 
                                                title = f'{columns_to_names[col]}',
                                                axs = [axs[c]])
        axs[c].grid(True)
    for c, ax in enumerate(axs):
        if c == 0:
            ax.set_ylabel("CDF")
        else:
            ax.set_ylabel("")
        
    fig.savefig('./plots/paper/fig5/change_in_rna_at_diff_megaloop_bins.pdf', bbox_inches='tight')



from volcano_plot import label_points
def plot_metadomain_score_fc_fc_plot(pval_df, stat_df, ind_to_gene):
    sign_cluster0 = (pval_df[0] < .01) * np.sign(stat_df[0]) * (stat_df[0].abs() > 1.5)
    sign_cluster4 = (pval_df[4] < .01) * np.sign(stat_df[4]) * (stat_df[4].abs() > 1.5)

    color = np.zeros_like(sign_cluster0).astype('str')
    color[(sign_cluster0 == 0) & (sign_cluster4 == 0)] = 'lightgray'
    color[(sign_cluster0 == 1) & (sign_cluster4 == 0)] = 'salmon'
    color[(sign_cluster0 == 0) & (sign_cluster4 == 1)] = 'pink'
    color[(sign_cluster0 == 1) & (sign_cluster4 == 1)] = 'red'
    color[(sign_cluster0 == -1) & (sign_cluster4 == 0)] = 'aqua'
    color[(sign_cluster0 == 0) & (sign_cluster4 == -1)] = 'skyblue'
    color[(sign_cluster0 == -1) & (sign_cluster4 == -1)] = 'steelblue'
    color[(sign_cluster0 == -1) & (sign_cluster4 == 1)] = 'fuchsia'
    color[(sign_cluster0 == 1) & (sign_cluster4 == -1)] = 'fuchsia'

    name_dict = {
        'lightgray' : 'NS both',
        'salmon' : 'Treg-up Constit.',
        'pink' : 'Treg-up Dynamic',
        'red' : 'Treg-up Both',
        'aqua' : 'Treg-down Constit.',
        'skyblue' : 'Treg-down Dynamic',
        'steelblue' : 'Treg-down Both',
        'fuchsia' : 'Discordant',
    }

    fig, axs = init_subplots_exact(1, 1, fgsz=(25*mm, 25*mm), dpi = 200)
    scatter_with_pearson(stat_df[4], stat_df[0], plt.gca(), s = 1, c = color)
    for color_name, name in name_dict.items():
        n_with_color = (color==color_name).sum()
        label = f'{name}: n = {n_with_color}'
        plt.scatter([], [], color=color_name, label=label, s = 4, linewidth=.1, edgecolor='black')

    inds_to_label = [278, 6736, 5699, 242, 6400, 5529, 3207]
    names = [get_name(x, ind_to_gene) if x!=242 else 'Ctla4' for x in inds_to_label]
    label_points(stat_df[4], stat_df[0], inds_to_label, names, s = 3, linewidth=.38, bold=False)

    plt.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False)

    axs.set_axisbelow(True)
    add_xaxis_labels('Tcon', 'Treg', axs, fontsize=6)
    add_yaxis_labels('Tcon', 'Treg', axs, fontsize=6)

    plt.xlim([-11, 11])
    plt.ylim([-11, 11])
    z = np.linspace(-11, 11)
    plt.plot(z, z, color='black', linestyle='--')
    axs.set_xticks([-5, 0, 5])
    axs.set_yticks([-5, 0, 5])
    plt.xlabel("Active 2")
    plt.ylabel("Active 1")
    plt.title("Change in O/E")
    fig.savefig('./plots/paper/fig5/change_in_oe_plot.pdf', bbox_inches='tight')




def get_weight(tmp_mat, u, v, node_to_ind):
    weight1, weight2 = tmp_mat[node_to_ind[u]].sum(), tmp_mat[node_to_ind[v]].sum()
    return weight1, weight2
def prune_tcon_edges(G, tcon_tmp, node_to_ind, cutoff=60):
    edges_to_remove = []
    for u, v in G.edges():
        weight = max(tcon_tmp[node_to_ind[u]].sum(), tcon_tmp[node_to_ind[v]].sum())
        if weight < cutoff:
            edges_to_remove.append((u, v))
    G.remove_edges_from(edges_to_remove)

def prune_treg_edges(G, treg_tmp, node_to_ind):
    edges_to_remove = []
    for u, v in G.edges():
        weight = max(treg_tmp[node_to_ind[u]].sum(), treg_tmp[node_to_ind[v]].sum())
        if weight < 25:
            edges_to_remove.append((u, v))
    G.remove_edges_from(edges_to_remove)

def rotate_pos(pos, angle_deg):
    """Rotate positions of nodes by a specified angle in degrees."""
    angle_rad = np.radians(angle_deg)
    rotation_matrix = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                                [np.sin(angle_rad),  np.cos(angle_rad)]])
    
    pos_rotated = {}
    for node, xy in pos.items():
        pos_rotated[node] = np.dot(rotation_matrix, xy)
    return pos_rotated


def scale_pos_y(pos, scale_y):
    """Scale positions of nodes in the Y direction by a specified factor."""
    pos_scaled = {}
    for node, (x, y) in pos.items():
        pos_scaled[node] = (x, y * scale_y)
    return pos_scaled

def scale_pos_x(pos, scale_x):
    """Scale positions of nodes in the Y direction by a specified factor."""
    pos_scaled = {}
    for node, (x, y) in pos.items():
        pos_scaled[node] = (x* scale_x, y)
    return pos_scaled

from matplotlib.colors import LinearSegmentedColormap
def networkx_cluster2_embedding(
        self, goodinds, ind_to_gene):
    goodinds = self.goodinds[np.isin(self.merged_clustdict['all'], [4])]

    # goodinds = goodinds[inter_and_intra_connections[goodinds, :][:, goodinds].sum(axis=1)>15]
    focal_submat = self.inter_and_intra_connections_treg[goodinds, :][:, goodinds].copy().astype(float)
    focal_submat2 = self.inter_and_intra_connections_tcon[goodinds, :][:, goodinds].copy().astype(float)
    final_focal_submat = np.zeros_like(focal_submat)

    final_focal_submat[(focal_submat == 0) & (focal_submat2 == 0)] = 0
    final_focal_submat[(focal_submat == 1) & (focal_submat2 == 1)] = 0
    final_focal_submat[(focal_submat == 1) & (focal_submat2 == 0)] = 1
    final_focal_submat[(focal_submat == 0) & (focal_submat2 == 1)] = -1




    # treg_submat = inter_oe_mat_treg[goodinds, :][:, goodinds].copy()
    # tcon_submat = inter_oe_mat_tcon[goodinds, :][:, goodinds].copy()
    delta_submat = np.concatenate([
                                # treg_submat - tcon_submat, 
                                # treg_submat, tcon_submat,
                                focal_submat, focal_submat2,
                                # focal_submat - focal_submat2,
                                ], axis=1)

    # delta_submat = inter_oe_mat_tcon[goodinds, :][:, goodinds].copy().astype(float)

    genesoi = [243, 278, 8298, 2178, 2179, 3207, 3208, 6155, 6154, 6501, 6502, 9162, 5217, 2281, 826, 3784,
            2281, 826, 3784, 2187, 1974, 2035, 2556, 2281, 826, 3784, 2187, 1974, 2035, 2556, 6891, 8327, 9522, 6892, 3511, 8047
            ] + list(goodinds[final_focal_submat.sum(axis=1) > 10]) + list(goodinds[final_focal_submat.sum(axis=1) < -40])
    genesoi = list(set(genesoi))
    genesoi.remove(4712)
    genesoi.remove(5074)
    genesoi.remove(3037)
    genesoi.remove(1721)
    genesoi.remove(1915)
    genesoi.remove(8205)
    genesoi.remove(9551)
    genesoi.remove(1974)
    genesoi.remove(3511)
    genesoi.remove(6891)
    genesoi.remove(8047)
    genesoi.remove(8327)
    genesoi.remove(9522)
    genesoi.remove(8298)
    genesoi.remove(4894)
    genesoi.remove(3606)
    genesoi.remove(2187)
    genesoi.remove(826)
    genesoi.remove(354)
    new_names = {}
    nodes_to_label = []
    n = len(goodinds)
    for i in np.arange(n):
        if goodinds[i] in genesoi:
            name = get_name(goodinds[i], ind_to_gene, default=goodinds[i])
            if name in nodes_to_label:
                name = name + "_2"
            if name == 'Bach2it1':
                name = 'Bach2'
            new_names[i] = name
            nodes_to_label.append(name)
        else:
            new_names[i] = goodinds[i]


    gist_heat_r = plt.cm.get_cmap('gist_heat_r', 256)
    gist_heat_r_colors = gist_heat_r(np.linspace(0, 1, 256))

    # Modifying the colors: replace red with blue
    blue_heat_colors = np.copy(gist_heat_r_colors)
    blue_heat_colors[:, 0] = gist_heat_r_colors[:, 2]  # Swap red channel with blue
    blue_heat_colors[:, 2] = gist_heat_r_colors[:, 0]  # Swap blue channel with red

    # Creating a new colormap
    blue_heat = LinearSegmentedColormap.from_list("blue_heat", blue_heat_colors)


    delta_submat = delta_submat.clip(-5, 5)
    delta_submat[np.isnan(delta_submat)] = 0
    # m = -delta_submat
    m = np.corrcoef(delta_submat)
    m = m - np.diag(np.diag(m))
    # co = .5
    co = .00001
    m[m < -co] = -co
    uco = .5
    m[m > uco] = uco

    G = nx.from_numpy_array(m)
    G = nx.relabel_nodes(G, new_names)
    n = len(final_focal_submat)
    np.random.seed(0)
    pos = nx.spring_layout(G)
    pos = rotate_pos(pos, -40)
    # Apply the scaling to the original positions
    pos = scale_pos_y(pos, 0.45)  # Scale Y by 0.5
    pos = scale_pos_x(pos, 0.75)  # Scale Y by 0.5

    # pos = nx.nx_agraph.graphviz_layout(G, prog="fdp")


    # pos = nx.kamada_kawai_layout(G)

    # pos = nx.force_directed(G)

    tcon_tmp = focal_submat2.copy()
    tcon_tmp[(focal_submat2==1) & (focal_submat==1)] = 0
    G = nx.from_numpy_array(tcon_tmp)
    G = nx.relabel_nodes(G, new_names)

    treg_tmp = focal_submat.copy()
    treg_tmp[(focal_submat==1) & (focal_submat2==1)] = 0

    assert len(goodinds) == len(G.nodes())

    high_weight_threshold = 0

    # Categorize edges and assign colors
    edge_colors = ['red' if G[u][v]['weight'] > high_weight_threshold else 'gray' for u, v in G.edges()]
    edge_colors = ['gray' if G[u][v]['weight'] > high_weight_threshold else 'gray' for u, v in G.edges()]
    node_to_ind = transpose_dict(dict(enumerate(G.nodes())))
    edge_colors = [blue_heat(tcon_tmp[node_to_ind[u]].sum()/200) if G[u][v]['weight'] > high_weight_threshold else 'gray' for u, v in G.edges()]

    # Draw the graph
    fig, axs = init_subplots_exact(1, 1, fgsz=(60*mm, 60*mm), dpi = 300)

    node_sizes = [tcon_tmp[c].sum() if node in nodes_to_label else 0 for c, node in enumerate(G.nodes())]
    node_sizes = arr([max(treg_tmp[c].sum()+1, tcon_tmp[c].sum()+1) for c, node in enumerate(G.nodes())])
    # node_sizes = arr([((treg_tmp[c] > 0) | (tmp[c]>0)).sum() for c, node in enumerate(G.nodes())])
    node_sizes[node_sizes < 2] = 0
    node_sizes /= 5
    node_sizes *= node_sizes/20

    node_colors = [
            sns.color_palette('coolwarm', as_cmap=True)(final_focal_submat[c].sum()/100+.5) for c, node in enumerate(G.nodes())]


    prune_tcon_edges(G, tcon_tmp, node_to_ind, cutoff=25)
    nx.draw(G, pos, with_labels=False, node_color=node_colors, edge_color=edge_colors, width=.02,
            node_size = node_sizes, font_size=2, linewidths=0, 
            #node_line='black',
            ax = axs
        ) 

    label_pos = {node: (position[0], position[1] + 0.02) for node, position in pos.items()}  # Slightly shift labels

    # Label only selected nodes
    idx = np.isin(list(G.nodes()), nodes_to_label)
    labels = {node: node for node in nodes_to_label}
    nx.draw_networkx_labels(G, label_pos, labels=labels, font_size=4)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, edgecolors='black', linewidths=0.1*0, node_size=node_sizes*1.25, ax=axs)
    nx.draw_networkx_nodes(G, pos, nodelist=nodes_to_label, 
                            node_color = 'None',
                        edgecolors='black', linewidths=0.2, node_size=arr(node_sizes)[idx]*1.25, ax=axs)


    G = nx.from_numpy_array(treg_tmp)
    G = nx.relabel_nodes(G, new_names)

    # pos = nx.graphviz_layout(G)
    # pos = nx.kamada_kawai_layout(G)

    high_weight_threshold = 0


    def get_weight_of_edge_treg(u, v):
        weight1, weight2 = get_weight(treg_tmp, u, v, node_to_ind)
        return (int(weight1 > 25) + int(weight2 > 25))/2

    # def get_weight_of_edge_treg2(u, v):
    #     weight1, weight2 = get_weight(treg_tmp, u, v)
    #     return weight1+weight2
    #     # return (max(weight1, weight2))

    # Categorize edges and assign colors
    edge_colors = ['red' if G[u][v]['weight'] > high_weight_threshold else 'gray' for u, v in G.edges()]
    edge_colors = ['gray' if G[u][v]['weight'] > high_weight_threshold else 'gray' for u, v in G.edges()]
    prune_treg_edges(G, treg_tmp, node_to_ind)
    # edge_colors = [sns.color_palette('gist_heat', as_cmap=True)(get_weight_of_edge(u, v)) if G[u][v]['weight'] > high_weight_threshold else 'gray' for u, v in G.edges()]
    # edge_colors = [sns.color_palette('gist_heat_r', as_cmap=True)(get_weight_of_edge_treg(u, v)*.25) for u, v in G.edges()]
    edge_colors = [sns.color_palette('gist_heat_r', as_cmap=True)(get_weight_of_edge_treg(u, v)*.25) for u, v in G.edges()]
    # edge_colors = [sns.color_palette('gist_heat_r', as_cmap=True)(get_weight_of_edge_treg2(u, v)/400) for u, v in G.edges()]

    nx.draw(G, pos, with_labels=False, node_color=node_colors, edge_color=edge_colors, width=.1,
            node_size=0, 
            ax = axs
        ) 
    plt.ylim([-.5, .5])
    plt.xlim([-1, 1])
    fig.savefig('./plots/paper/fig5//interchrom_hub_as_network.pdf', bbox_inches='tight')    


from plotting_functions import add_pval_to_plot
from aux_functions import nonan_test


def hub_atac_coaccessibility():
    atac_coac_pref = './coaccessibility_analysis/'
    # all_cluster_atac = pd.read_parquet(atac_coac_pref + 'all_treg_interchromosomal_metadomain_cluster_corr.parquet')
    # refined_metadomain_df = pd.read_csv('final_loops/metaloops/refined_metaloops/refined_metaloop_anchors.bed', sep='\t', header=None)

    tmpdf1_intra = pd.read_csv('intermediate_files/intrachromosomal_treg_metadomain_bedfile.bed', sep='\t', header=None)
    anc1 = tmpdf1_intra[[0, 1, 2]]
    anc2 = tmpdf1_intra[[3, 4, 5]]
    anc2.columns = anc1.columns
    # intra_metadomain_anchors = pd.concat([anc1, anc2], axis=0).value_counts()

    tmpdf1_inter = pd.read_csv('intermediate_files/all_treg_interchromosomal_metadomain_bedfile.bed', sep='\t', header=None)
    anc1 = tmpdf1_inter[[0, 1, 2]]
    anc2 = tmpdf1_inter[[3, 4, 5]]
    anc2.columns = anc1.columns
    # inter_metadomain_anchors = pd.concat([anc1, anc2], axis=0).value_counts()

    # intra_metadomain_anchors_bedtool = pbt.BedTool.from_dataframe(intra_metadomain_anchors.reset_index())
    # inter_metadomain_anchors_bedtool = pbt.BedTool.from_dataframe(inter_metadomain_anchors.reset_index())

    # _5kb_anchors_bedtool = pbt.BedTool.from_dataframe(refined_metadomain_df.value_counts().reset_index())
    # all_50kb_metadomain_anchors = intra_metadomain_anchors_bedtool.cat(inter_metadomain_anchors_bedtool)

    all_5kb_bins_atac_coac = pd.read_parquet(atac_coac_pref + 'all_5kb_anchors_in_50kb_Treg_metadomains_atac_corr.parquet')

    vs = np.triu(all_5kb_bins_atac_coac, k=1)
    vs[vs==0] = np.nan
    vs[~np.isnan(vs)] = 0
    all_5kb_bins_atac_coac += vs
    def get_anchor_cluster_annotation(coldata):
        rows = []
        og_vals = []
        for x in coldata.unique():
            chrom, s = x.split(':')
            rows.append([chrom, int(s), int(s) + 5000])
            og_vals.append(x)
        col_rows = pbt.BedTool(rows)
        col_rows_with_annot = col_rows.intersect(add_chr_to_bedtool(pbt.BedTool('for_susie/250kb_anchors_with_cluster_label.bed')), wao=True)
        cluster_annotation = get_col(col_rows_with_annot, -2)
        cluster_annotation = cluster_annotation.astype('<U16')
        cluster_annotation[cluster_annotation == '.'] = '-1'
        cluster_annotation = cluster_annotation.astype(int)
        cluster_annotation_dict = dict(zip(og_vals, cluster_annotation))
        print("Done making dict!")
        cluster_annotation = [cluster_annotation_dict[x] for x in coldata]
        return cluster_annotation

    # def get_anchor_metadomain_annotation(coldata):
    #     rows = []
    #     og_vals = []
    #     for x in coldata.unique():
    #         chrom, s = x.split(':')
    #         rows.append([chrom, int(s), int(s) + 5000])
    #         og_vals.append(x)
    #     col_rows = pbt.BedTool(rows)
    #     col_rows_with_annot = col_rows.intersect(add_chr_to_bedtool(pbt.BedTool('for_susie/250kb_anchors_with_cluster_label.bed')), wao=True)
    #     cluster_annotation = get_col(col_rows_with_annot, -2)
    #     cluster_annotation = cluster_annotation.astype('<U16')
    #     cluster_annotation[cluster_annotation == '.'] = '-1'
    #     cluster_annotation = cluster_annotation.astype(int)
    #     cluster_annotation_dict = dict(zip(og_vals, cluster_annotation))
    #     print("Done making dict!")
    #     cluster_annotation = [cluster_annotation_dict[x] for x in coldata]
    #     return cluster_annotation

    def susie_anchor_to_grange(susie_anchor):
        chrom, s = susie_anchor.split(":")
        s = int(s)
        s = s//50_000*50_000
        e = s + 50_000
        anc1 = tuple_to_grange(chrom, s, e)
        return anc1


    melted_atac_coac = all_5kb_bins_atac_coac.iloc[::1, ::1].unstack().reset_index()
    melted_atac_coac.columns = 'anchor1', 'anchor2', 'r'
    melted_atac_coac = melted_atac_coac.dropna()

    melted_atac_coac['cluster1'] = get_anchor_cluster_annotation(melted_atac_coac['anchor1'])
    melted_atac_coac['cluster2'] = get_anchor_cluster_annotation(melted_atac_coac['anchor2'])
    subdf = melted_atac_coac[(melted_atac_coac['cluster1'].isin([0, 2, 9])) 
                        & (melted_atac_coac['cluster2'].isin([0, 2, 9]))]

    loop_dict = {}
    for _, i in pd.concat([tmpdf1_intra, tmpdf1_inter], axis=0).iterrows():
        anc1, anc2 = i[:3], i[3:]
        grange1 = tuple_to_grange(*add_chr_to_anc(tuple(anc1)))
        grange2 = tuple_to_grange(*add_chr_to_anc(tuple(anc2)))
        loop_dict[(grange1, grange2)] = 1
        loop_dict[(grange2, grange1)] = 1

    import warnings

    warnings.filterwarnings(category=FutureWarning, action='ignore')

    has_metadomain = []
    for i1, i2 in zip(subdf['anchor1'], subdf['anchor2']):
        anc1, anc2 = susie_anchor_to_grange(i1), susie_anchor_to_grange(i2)
        if (anc1, anc2) in loop_dict:
            has_metadomain.append(True)
        else:
            has_metadomain.append(False)

    subdf['has_metadomain'] = has_metadomain

    colors = ['lightgreen', 'green', 'red']
    xlims = [[-.2, .2], [-.4, .6]]
    titles = ['Active 1', 'Active 2']
    fig, axs = init_subplots_exact(2, 1, fgsz=(30*mm, 30*mm), dpi = 200, sharey=True)
    for c, (u1, u2) in enumerate([[0, 0], [2, 2]]):
        same_cluster = ((subdf['cluster1'] == u1) & (subdf['cluster2'] == u2)) | ((subdf['cluster1'] == u2) & (subdf['cluster2'] == u1))
        has_metadomain = subdf['has_metadomain'] > 0
        idx_mega = has_metadomain & same_cluster
        idx_nonmega = (~has_metadomain) & same_cluster
        plt.sca(axs[c])
        print(scipy.stats.ranksums(subdf.loc[idx_mega, 'r'], subdf.loc[idx_nonmega, 'r']))
        sns.ecdfplot(subdf.loc[idx_mega, 'r'], label=f'With metadomain', c=colors[c])
        sns.ecdfplot(subdf.loc[idx_nonmega, 'r'], label=f'Without metadomain', c=colors[c], linestyle='--')
        plt.xlim(xlims[c])
        if c > 0:
            plt.ylabel("")
        plt.xlabel("Pearson Correlation")
        plt.title(titles[c])
        plt.legend(bbox_to_anchor=(0, .95), loc='upper left', frameon=False,
                handlelength=1.2)
        add_pval_to_plot(subdf.loc[idx_mega, 'r'], subdf.loc[idx_nonmega, 'r'], plt.gca(),
                        test=nonan_test)
    fig.savefig('./plots/paper/fig6/atac_hub_coacc.pdf', bbox_inches='tight')


def metaloop_coaccessibility():
    df1 = pd.read_csv('coaccessibility_analysis/metaloop_corr_anchor_v_anchor.txt', sep='\t')
    df2 = pd.read_csv('coaccessibility_analysis/metaloop_corr_anchor_v_anchor_random.txt', sep='\t')
    df3 = pd.read_csv('coaccessibility_analysis/metaloop_corr_nonanchor_random.txt', sep='\t')

    fig, axs = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi = 200)
    sns.ecdfplot(df1['r'], color='purple',  label='Metaloop')
    sns.ecdfplot(df2['r'], color='purple', linestyle='--', label='Random pairs \nfrom same anchors')
    plt.xlabel("Correlation")
    plt.title("ATAC correlation, metaloops")
    plt.legend(frameon=False, loc='lower right', handlelength=1,  fontsize=6)
    plt.xlim([-.2, .3])
    add_pval_to_plot(df1['r'], df2['r'], plt.gca(), test=nonan_test)
    fig.savefig('./plots/paper/fig6/atac_metaloop_coacc.pdf', bbox_inches='tight')    


def make_h3k27ac_stat5_motif_plot(h3k27ac_motif_df, motif='Stat5a'):
    all_h3k27ac = add_chr_to_bedtool(pbt.BedTool('peaks/differential/all_threshold_27ac.csv')).to_dataframe()
    all_h3k27ac['score'] = -all_h3k27ac['score']
    all_h3k27ac = all_h3k27ac[(all_h3k27ac['name'] < 1000) & (all_h3k27ac['name'] > 10)]
    all_h3k27ac.index = bedtool_to_index(pbt.BedTool.from_dataframe(all_h3k27ac))

    outpref = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/rudensky_scrna/prelim-analysis/call_motifs/bedfiles/h3k27ac_treg.csv'
    df = add_chr_to_bedtool(pbt.BedTool('peaks/differential/all_threshold_27ac.csv')).to_dataframe().iloc[:, :3]
    pbt.BedTool.from_dataframe(df).saveas(outpref)

    fig, axs = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi = 200)

    v_dict = {}

    h3k27ac_motif_df = h3k27ac_motif_df.loc[all_h3k27ac.index]
    n_motifs = h3k27ac_motif_df[motif].clip(0, 2)

    kwargs = [dict(color='gray', linestyle='--'), dict(color='salmon'), dict(color='red'), dict(color='maroon')]
    for c, u in enumerate(np.unique(n_motifs)):
        idx = (n_motifs == u).values
        v = all_h3k27ac.loc[idx, 'score']
        if u == 0:
            sns.ecdfplot(v, **kwargs[c], label=f'{u} Stat5 Motifs')
        else:
            p = format_pvalue(scipy.stats.ranksums(v, v_dict[0])[1])
            sns.ecdfplot(v, **kwargs[c], label=f'{u} Stat5 Motifs;\np={p}')
        v_dict[u] = v
        
    plt.xlim([-2, 2])
    plt.legend(loc='upper left', frameon=False)
    add_xaxis_labels("Tcon", "Treg", axs, fontsize=6)
    plt.xlabel("H3K27ac LFC")
    plt.title("H3K27ac LFC at Stat5 sites")
    plt.legend(bbox_to_anchor=(0, -.25), loc='upper left', ncol=3)
    plt.ylabel("CDF")
    fig.savefig('./plots/paper/fig7/stat5_motifs.pdf', bbox_inches='tight')    

def make_h3k27ac_stat5_motif_plot_1(motif='Stat5a'):
    chromvar_pref = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/rudensky_scrna/prelim-analysis/call_motifs/output/counts/h3k27ac_treg_Mus_musculus//h3k27ac_treg_p=1e-05.csv'
    chromvar_motifs = pd.read_csv(chromvar_pref, index_col = 0)
    all_h3k27ac = add_chr_to_bedtool(pbt.BedTool('peaks/differential/all_threshold_27ac.csv')).to_dataframe()
    all_h3k27ac['score'] = -all_h3k27ac['score']
    pick_idx = (all_h3k27ac['name'] < 1000) & (all_h3k27ac['name'] > 10)
    all_h3k27ac = all_h3k27ac[pick_idx]
    chromvar_motifs = chromvar_motifs.loc[pick_idx.values]
    
    outpref = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/rudensky_scrna/prelim-analysis/call_motifs/bedfiles/h3k27ac_treg.csv'
    df = add_chr_to_bedtool(pbt.BedTool('peaks/differential/all_threshold_27ac.csv')).to_dataframe().iloc[:, :3]
    pbt.BedTool.from_dataframe(df).saveas(outpref)

    fig, axs = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi = 200)
    
    v_dict = {}
    n_motifs = chromvar_motifs[motif].clip(0, 2)
    
    kwargs = [dict(color='gray', linestyle='--'), dict(color='salmon'), dict(color='red'), dict(color='maroon')]
    for c, u in enumerate(np.unique(n_motifs)):
        idx = (n_motifs == u).values
        v = all_h3k27ac.loc[idx, 'score']
        if u == 0:
            sns.ecdfplot(v, **kwargs[c], label=f'{u} Stat5 Motifs')
        else:
            p = format_pvalue(scipy.stats.ranksums(v, v_dict[0])[1])
            sns.ecdfplot(v, **kwargs[c], label=f'{u} Stat5 Motifs;\np={p}')
        v_dict[u] = v
        
    plt.xlim([-2, 2])
    plt.legend(loc='upper left', frameon=False)
    add_xaxis_labels("Tcon", "Treg", axs, fontsize=6)
    plt.xlabel("H3K27ac LFC")
    plt.title("H3K27ac LFC at Stat5 sites")
    plt.legend(bbox_to_anchor=(0, -.25), loc='upper left', ncol=3)
    fig.savefig('./plots/paper/fig7/stat5_motifs.pdf', bbox_inches='tight')    


def generate_all_bw_values_chrom(bw, chromsizes, parsed_chroms):
    all_vals = []
    for chrom in parsed_chroms:
        e = chromsizes[f'{chrom}']
        vals = bw.fetch(f'chr{chrom}', 0, e, bins = e//5_000)
        all_vals += list(np.ravel(vals))
    return arr(all_vals)

import bbi
def make_stat5_h3k27ac_dataframes(bbi_path_dict, chromsizes, parsed_chroms):
    keys = [
    ['Treg Stat5', 'Tcon Stat5'],
    ['Treg Ets1', 'Tconv Ets1'],
    ['Treg CREB', 'Tconv CREB'],
    ['Treg Satb1', 'Tconv Satb1'],
    ['Treg CTCF', 'CD4 CTCF'],
    ['Treg Smc1a', 'Tcon Smc1a'],
    ['Treg H3K27ac', 'Tcon H3K27ac'],
    ['Treg Tet2', 'Tcon Tet2'],

    ['Treg Med1', 'Tcon Med1'], 
    ['Treg Runx1', 'Tcon Runx1'], 
    ['Treg Bcl11b', 'Tcon Bcl11b'], 
    ['Treg Ep300', 'Tcon Ep300'],
    ['Treg Foxp3', None]
    ]

    baseline_df = pd.DataFrame()
    val_df = pd.DataFrame()
    for key1, key2 in keys:
        name = key1.split(" ")[1]
        if key1 == 'Treg Foxp3':
            v_treg = np.log2(1 + generate_all_bw_values_chrom(bbi.open(bbi_path_dict[key1]), chromsizes, parsed_chroms))
            tmp_delta_stat = (v_treg)
            val_df[name] = tmp_delta_stat
            baseline_df[name] = v_treg
        else:
            v_treg = np.log2(1 + generate_all_bw_values_chrom(bbi.open(bbi_path_dict[key1]), chromsizes, parsed_chroms))
            v_tcon = np.log2(1 + generate_all_bw_values_chrom(bbi.open(bbi_path_dict[key2]), chromsizes, parsed_chroms))
            tmp_delta_stat = (v_treg - v_tcon)
            name = key1.split(" ")[1]
            val_df[name] = tmp_delta_stat
            baseline_df[name] = v_treg    
    data = val_df[~(val_df == 0).all(axis=1)]
    data = data[['H3K27ac', 'Stat5', 'Smc1a', 'Tet2', 'Ets1', 'Bcl11b', 'Runx1', 'Satb1', 'CREB', 'CTCF', 'Foxp3']]

    baseline_data = baseline_df[~(baseline_df == 0).all(axis=1)]
    baseline_data = baseline_data[['H3K27ac', 'Stat5', 'Smc1a', 'Tet2', 'Ets1', 'Bcl11b', 'Runx1', 'Satb1', 'CREB', 'CTCF', 'Foxp3']]
    return baseline_data, data

def stat5_h3k27ac_barplots(baseline_data, data):
    order = ['Stat5', 'Smc1a', 'Tet2', 'Ets1', 'Bcl11b', 'Runx1',  'Foxp3', 'Satb1', 'CREB', 'CTCF']
    fig, axs = init_subplots_exact(1, 1, fgsz=(10*mm, 30*mm), dpi = 200, sharex=True, as_list=True)
    plt.sca(axs[0])
    tmp = baseline_data.corr().loc['H3K27ac'].drop("H3K27ac").loc[order]
    tmp.plot.barh(zorder=3, width=.75)
    plt.xticks(rotation=0)
    plt.ylabel("")
    plt.title("Correlation with H3K27ac")
    plt.xlim([0, 1])
    plt.xlabel("Correlation")
    plt.gca().grid(False)
    fig.savefig('./plots/paper/fig7/h3k27ac_base_corr.pdf', bbox_inches='tight')


    order = ['Stat5', 'Smc1a', 'Tet2', 'Ets1', 'Bcl11b', 'Runx1', 'Foxp3', 'Satb1', 'CREB', 'CTCF']
    fig, axs = init_subplots_exact(1, 1, fgsz=(10*mm, 30*mm), dpi = 200, sharex=True, as_list=True)
    plt.sca(axs[0])
    tmp = data.corr().loc['H3K27ac'].drop("H3K27ac").loc[order]
    tmp.plot.barh(zorder=3, width=.75)
    plt.xticks(rotation=0)
    plt.ylabel("")
    plt.title("Correlation with H3K27ac")
    plt.xlim([-.2, 1])
    plt.xlabel("Correlation")
    plt.gca().grid(False)
    fig.savefig('./plots/paper/fig7/h3k27ac_lfc_corr.pdf', bbox_inches='tight')    


def stat5_h3k27ac_scatterplot(baseline_data):
    fig, axs = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi = 200)
    plt.scatter(baseline_data['Stat5'], baseline_data['H3K27ac'],  s=.01, rasterized=True)
    add_pearson_to_plot(baseline_data['Stat5'], baseline_data['H3K27ac'], plt.gca(), test=samesize_nonan_test)
    plt.xlabel("Stat5")
    plt.ylabel("H3K27ac")
    plt.title("Stat5 - H3K27ac Correlation")
    axs.set_axisbelow(True)
    fig.savefig('./plots/paper/fig7/stat_h3k27ac_corr.pdf', bbox_inches='tight', dpi=2000)


def chip_hub_correlation_plots(self, bw_val_df_all_250kb, stat_df, columns_to_names, row_colors_dict):
    delta_stat = {
        'H3K27ac' : np.log2(np.ravel(bw_val_df_all_250kb['Treg H3K27ac'] / bw_val_df_all_250kb['Tcon H3K27ac'])),
        'Stat5' : np.log2(np.ravel(bw_val_df_all_250kb['Treg Stat5'] / bw_val_df_all_250kb['Tcon Stat5'])),
        'Smc1a' : np.log2(np.ravel(bw_val_df_all_250kb['Treg Smc1a'] / bw_val_df_all_250kb['Tcon Smc1a'])),
        'Ets1' : np.log2(np.ravel(bw_val_df_all_250kb['Treg Ets1'] / bw_val_df_all_250kb['Tconv Ets1'])),
        'CTCF' : np.log2(np.ravel(bw_val_df_all_250kb['Treg CTCF'] / bw_val_df_all_250kb['CD4 CTCF'])),
        'CREB' : np.log2(np.ravel(bw_val_df_all_250kb['Treg CREB'] / bw_val_df_all_250kb['Tconv CREB'])),
        'Satb1' : np.log2(np.ravel(bw_val_df_all_250kb['Treg Satb1'] / bw_val_df_all_250kb['Tconv Satb1'])),
        'Foxp3' : np.log2(np.ravel(bw_val_df_all_250kb['Treg Foxp3'])),
        'H3K4me1' : np.log2(np.ravel(bw_val_df_all_250kb['Treg H3K4me1'] / bw_val_df_all_250kb['Tcon H3K4me1'])),
        'H3K4me3' : np.log2(np.ravel(bw_val_df_all_250kb['Treg H3K4me3'] / bw_val_df_all_250kb['Tcon H3K4me3'])),
        'H3K27me3' : np.log2(np.ravel(bw_val_df_all_250kb['Treg H3K27me3'] / bw_val_df_all_250kb['Tcon H3K27me3'])),
    }

    row_colors = ['lightgreen', 'green', 'orange']

    r_df = []
    for name, delta_vector in delta_stat.items():
        rs = []
        inds = []
        for cluster in self.merged_inds_to_subset:
            indsoi = self.goodinds[self.merged_clustdict['all'] == cluster]
            if indsoi.sum() == 0:
                continue
            subdf = stat_df.loc[:, cluster]
            x, y = subdf, delta_vector[subdf.index]
            r = samesize_nonan_test(x, y)[0]
            rs.append(r)
            inds.append(cluster)
        r_df.append(rs)
    r_df = pd.DataFrame(r_df, index=delta_stat.keys())
    r_df.columns = [columns_to_names[x] for x in inds]
    print(2)
    gc = sns.clustermap(r_df, col_cluster=False, row_cluster=False, 
                        figsize=(58*mm, 58*mm), col_colors=row_colors,
                        annot=True, cmap='coolwarm', fmt='.2f', 
                        cbar_pos=(.9, .8, .05, .1), vmin=-.7, vmax=.7,
                        annot_kws = {'fontsize' : 4}, zorder=3,
                    )

    gc.fig.dpi = 300
    col_colors_pos = gc.ax_col_colors.get_position().bounds
    new_height = 0.05  # Adjust this value to your liking
    gc.ax_col_colors.set_position([col_colors_pos[0], col_colors_pos[1], col_colors_pos[2], new_height])
    gc.ax_heatmap.set_title('ChIP - Hub\nCorrelation (LFC)', y = 1.1)

        
    gc.ax_heatmap.set_xlabel("Cluster")
    gc.ax_heatmap.set_ylabel("Factor",)
    plt.sca(gc.ax_heatmap)
    ticks = []
    for i, ticklabel in enumerate(gc.ax_heatmap.get_xticklabels()):
        ticklabel.set_color(row_colors_dict[ticklabel.get_text()])
        t = ticklabel.get_text()
        t = t.replace(")", '').replace("(", '')
        ticks.append(t)
    plt.gca().set_xticklabels(ticks)
    plt.xticks(rotation=0, fontsize=6)

    gc.savefig('./plots/paper/fig7/chip_binding_to_metaloops.pdf', bbox_inches='tight')