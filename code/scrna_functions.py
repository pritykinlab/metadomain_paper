import numpy as np
arr = np.asarray
def get_distances(anclist):
    dists = arr(anclist)[:, 1].astype(int)
    return dists

def calculate_loop_corr_5kb(cormat, tr_adata, gene2anc, megaloop_subloops):
    megaloop_subloops = set(megaloop_subloops)
    arr = np.asarray
    has_megaloop = []
    no_megaloop = []
    has_megaloop_names = []
    no_megaloop_names = []
    for c1, gene1 in enumerate(tr_adata.var.index):
        for c2, gene2 in enumerate(tr_adata.var.index):
            if gene1 <= gene2:
                continue
            r = cormat[c1, c2]
            anc1, anc2 = gene2anc.get(gene1, None), gene2anc.get(gene2, None)
            if (anc1 is None) or (anc2 is None):
                continue
            dist1, dist2 = get_distances(anc1), get_distances(anc2)
            distances = np.subtract.outer(arr(dist1), arr(dist2))/1e6
            loopcount = 0
            for a1 in anc1:
                for a2 in anc2:
                    fullloop = tuple(list(a1) + list(a2))
                    transpose_fullloop = tuple(list(a2) + list(a1))
                    # print(fullloop, list(megaloop_subloops)[0])
                    # raise Exception
                    if (fullloop in megaloop_subloops) or (transpose_fullloop in megaloop_subloops):
                            loopcount += 1
            if loopcount > 0:
                has_megaloop.append(r)
                has_megaloop_names.append([gene1, gene2])
            elif loopcount == 0:
                no_megaloop.append(r)
                no_megaloop_names.append([gene1, gene2])
    return has_megaloop, no_megaloop, has_megaloop_names, no_megaloop_names
import time
# def calculate_megaloops_5kb(cormat, tr_adata, gene2adata_ind, anc2gene, gene2dist, megaloop_subloops, ):
#     megaloop_subloops = set(megaloop_subloops)
#     genes_have_loop = np.zeros_like(cormat)
#     for loop in megaloop_subloops:
#         loop1, loop2 = loop[:3], loop[3:6]
#         genes1, genes2 = anc2gene.get(loop1, []), anc2gene.get(loop2, [])
#         for g1 in genes1:
#             for g2 in genes2:
#                 i1 = gene2adata_ind.get(g1, None)
#                 i2 = gene2adata_ind.get(g2, None)
#                 if (i1 is None) or (i2 is None):
#                     continue
#                 i1, i2 = sorted([i1, i2])
#                 genes_have_loop[i1, i2] = 1
#                 genes_have_loop[i2, i1] = 1
#     good = np.where(np.triu(genes_have_loop, k=1))
#     bad = np.where(np.triu(1-genes_have_loop, k=1))
#     return genes_have_loop, cormat[good], cormat[bad], 

import scipy
import scipy.stats
def compute_megaloop_enrichment_in_high_correlation(genematdict, adatadict, loop_bedfile, gene2dist):
    print("N anchors:", len(loop_bedfile))
    for c, cond in enumerate(genematdict):
        genemat = genematdict[cond].copy()
        adata = adatadict[cond].copy()
        
        names = adata.var.index
        cormat = np.corrcoef(genemat)
        
        gene2adata_ind = get_gene2adata_ind(adata)
        gene2anc, anc2gene, all_ancs = get_dict_from_loops(loop_bedfile)

        genes_have_loop, has_loop, no_loop = calculate_megaloops_5kb_fast(cormat, adata, gene2adata_ind, anc2gene, gene2dist, loop_bedfile)
        v1, v2 = (genes_have_loop[(np.triu(cormat>.1, k=1))], 
                genes_have_loop[(np.triu((np.abs(cormat)<0.1), k=1))])
        
        print(cond, scipy.stats.fisher_exact([[v1.sum(), (1-v1).sum()], 
                                [v2.sum(), (1-v2).sum()]]))    

import matplotlib.pyplot as plt
import seaborn as sns
from tad_functions import make_df_from_dict
def compute_corr_enrichment_in_megaloops(genematdict, adatadict, loop_bedfile, gene2dist, adata_type='HVG'):
    resdict = {}
    fig, axs = plt.subplots(1, 3, figsize=(8, 2))
    for c, cond in enumerate(genematdict):
        ax = axs[c]
        genemat = genematdict[cond].copy()
        adata = adatadict[cond].copy()
        
        v = np.ravel((adata.layers['poo']==0).mean(axis=0))
        good = v > -1
        adata = adata[:, good]
        genemat = genemat[good, :]

        names = adata.var.index
        cormat = np.corrcoef(genemat)
        
        gene2adata_ind = get_gene2adata_ind(adata)
        gene2anc, anc2gene, all_ancs = get_dict_from_loops(loop_bedfile)

        _, has_loop, no_loop = calculate_megaloops_5kb_fast(cormat, adata, gene2adata_ind, anc2gene, gene2dist, loop_bedfile)
        print(cond, scipy.stats.ranksums(has_loop[~np.isnan(has_loop)], no_loop[::1000][~np.isnan(no_loop[::1000])]))
        d = {'Has Loop' : list(has_loop),
            'No Loop' : list(no_loop[::100]),
            }
        df = make_df_from_dict(d)
        sns.boxplot(data=df, x='labels', y='values', ax=ax, fliersize=0)
        ax.set_title(f"Pearson r for {adata_type} \n with loops ({cond})")
        ax.set_xlabel("Label")
        ax.set_ylabel("Pearson R")
        ax.set_ylim([-.1, .1])
        resdict[cond] = d
    return fig, axs, resdict


def compute_corr_enrichment_in_megaloops(genematdict, adatadict, loop_bedfile, gene2dist, adata_type='HVG'):
    resdict = {}
    fig, axs = plt.subplots(1, 3, figsize=(8, 2))
    for c, cond in enumerate(genematdict):
        ax = axs[c]
        genemat = genematdict[cond].copy()
        adata = adatadict[cond].copy()
        
        #v = np.ravel((adata.layers['poo']==0).mean(axis=0))
        #good = v > -1
        #adata = adata[:, good]
        #genemat = genemat[good, :]

        names = adata.var.index
        cormat = np.corrcoef(genemat)
        
        gene2adata_ind = get_gene2adata_ind(adata)
        gene2anc, anc2gene, all_ancs = get_dict_from_loops(loop_bedfile)

        _, has_loop, no_loop = calculate_megaloops_5kb_fast(cormat, adata, gene2adata_ind, anc2gene, gene2dist, loop_bedfile)
        print(cond, scipy.stats.ranksums(has_loop[~np.isnan(has_loop)], no_loop[::1000][~np.isnan(no_loop[::1000])]))
        d = {'Has Loop' : list(has_loop),
            'No Loop' : list(no_loop[::100]),
            }
        df = make_df_from_dict(d)
        sns.boxplot(data=df, x='labels', y='values', ax=ax, fliersize=0)
        ax.set_title(f"Pearson r for {adata_type} \n with loops ({cond})")
        ax.set_xlabel("Label")
        ax.set_ylabel("Pearson R")
        ax.set_ylim([-.1, .1])
        resdict[cond] = d
    return fig, axs, resdict



def compute_corr_enrichment_in_anchors(genematdict, adatadict, loop_bedfile, gene2dist, adata_type='HVG', L = 0, R = np.inf):
    resdict = {}
    fig, axs = plt.subplots(1, 3, figsize=(8, 2))
    for c, cond in enumerate(genematdict):
        ax = axs[c]
        genemat = genematdict[cond].copy()
        adata = adatadict[cond].copy()

        names = adata.var.index
        cormat = np.corrcoef(genemat)
        
        gene2adata_ind = get_gene2adata_ind(adata)
        gene2anc, anc2gene, gene2chrom, all_loops = get_dict_from_loops(loop_bedfile)
        all_ancs = get_unique(loops_to_anchors(all_loops))

        _, has_loop, no_loop = calculate_megaloops_from_anchors(cormat, adata, gene2adata_ind, anc2gene, 
                                                                gene2dist, gene2chrom, all_ancs, L = L, R = R)

        d = {'Has Loop' : list(has_loop),
            'No Loop' : list(no_loop[::100]),
            }
        df = make_df_from_dict(d)
        sns.boxplot(data=df, x='labels', y='values', ax=ax, fliersize=0)
        ax.set_title(f"Pearson r for {adata_type} \n with loops ({cond})")
        ax.set_xlabel("Label")
        ax.set_ylabel("Pearson R")
        ax.set_ylim([-.1, .1])
        resdict[cond] = d
    return fig, axs, resdict

import sklearn
import sklearn.preprocessing
from sklearn.preprocessing import LabelEncoder
def calculate_megaloops_from_anchors(cormat, tr_adata, gene2adata_ind, anc2gene, gene2dist, gene2chrom, anchors,
                                        L = 0, R = np.inf):
    n = cormat.shape[0]
    genes_have_anchor = np.zeros(n)
    for anc in anchors:
        genes1 = anc2gene.get(anc, [])
        for g1 in genes1:
            i1 = gene2adata_ind.get(g1, None)
            if (i1 is None):
                continue
            genes_have_anchor[i1] = 1

    gene_locations = np.zeros(n)
    gene_chroms = np.zeros(n).astype(str)
    for gene, adata_ind in gene2adata_ind.items():
        location = gene2dist.get(gene, None)
        chrom = gene2chrom.get(gene, None)
        gene_locations[adata_ind] = location
        gene_chroms[adata_ind] = chrom

    le = LabelEncoder()
    chromvals = arr(le.fit_transform(gene_chroms))

    pairwise_chroms = np.equal.outer(chromvals, chromvals)
    pairwise_dists = np.abs(np.subtract.outer(gene_locations, gene_locations))
    genes_both_with_anchor = np.outer(genes_have_anchor, genes_have_anchor) > 0
    inds_to_compare = (pairwise_dists > L) & (pairwise_dists < R) & (pairwise_chroms) & (genes_both_with_anchor)
    inds_to_compare = np.triu(inds_to_compare + inds_to_compare.T, k=1)
    return genes_have_anchor, cormat[inds_to_compare], cormat[~inds_to_compare]


def calculate_megaloops_5kb_fast(cormat, tr_adata, gene2adata_ind, anc2gene, gene2dist, megaloop_subloops, ):
    megaloop_subloops = set(megaloop_subloops)
    genes_have_loop = np.zeros_like(cormat)
    for loop in megaloop_subloops:
        loop1, loop2 = loop[:3], loop[3:6]
        genes1, genes2 = anc2gene.get(loop1, []), anc2gene.get(loop2, [])
        for g1 in genes1:
            for g2 in genes2:
                i1 = gene2adata_ind.get(g1, None)
                i2 = gene2adata_ind.get(g2, None)
                if (i1 is None) or (i2 is None):
                    continue
                i1, i2 = sorted([i1, i2])
                genes_have_loop[i1, i2] = 1
    good = np.where(np.triu(genes_have_loop, k=1))
    bad = np.where(np.triu(1-genes_have_loop, k=1))
    return genes_have_loop, cormat[good], cormat[bad]

import aux_functions
from aux_functions import *

def get_dict_from_loops(loops):
    all_loops = set(loops)
    anchors = [make_str(l) for l in loops_to_anchors(all_loops)]
    gene2anc = {}
    anc2gene = {}
    peak_genes = pbt.BedTool('../peaks/RNA_coverage.narrowPeak')
    for i in peak_genes.intersect(anchors, wo=True):
        name = i[3]
        anc = i[-4:-1]
        anc = tuple(anc)
        anc2gene.setdefault(anc, [])
        gene2anc.setdefault(name, [])
        gene2anc[name].append(make_str(anc))    
        anc2gene[anc].append(name)
    
    gene2chrom = {}
    peak_genes = pbt.BedTool('../peaks/RNA_coverage.narrowPeak')
    for x in peak_genes:
        name = x[3]
        gene2chrom[name] = x[0]    
    return gene2anc, anc2gene, gene2chrom, all_loops


def get_gene2adata_ind(adata):
    names = adata.var.index
    n = len(names)
    return dict(zip(names, np.arange(n)))

def make_gene_mat_dict(adata_batch2):
    genematdict = {}
    adatadict = {}
    cells = ((adata_batch2.obs['ClusterCCA_Revised_Annotation_General'] == 'Treg') & 
            (adata_batch2.obs['sample_type'] == 'WT.TR'))

    treg_adata = adata_batch2[cells, :].copy()
    genemat = treg_adata.X.T
    subgenemat = genemat.copy()
    genematdict['Treg'] = subgenemat
    adatadict['Treg'] = treg_adata

    cells = ((adata_batch2.obs['ClusterCCA_Revised_Annotation_General'] == 'Tconv_resting') & 
            (adata_batch2.obs['sample_type'] == 'WT.TC'))
    tcon_adata = adata_batch2[cells, :].copy()
    genemat = tcon_adata.X.T
    subgenemat = genemat.copy()
    genematdict['Tconv_resting'] = subgenemat
    adatadict['Tconv_resting'] = tcon_adata

    cells = ((adata_batch2.obs['ClusterCCA_Revised_Annotation_General'] == 'Tconv_activated') & 
            (adata_batch2.obs['sample_type'] == 'WT.TC'))
    tcon_adata = adata_batch2[cells, :].copy()
    genemat = tcon_adata.X.T
    subgenemat = genemat.copy()
    genematdict['Tconv_activated'] = subgenemat
    adatadict['Tconv_activated'] = tcon_adata
    return genematdict, adatadict