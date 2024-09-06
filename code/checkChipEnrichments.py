import matplotlib.lines as mlines
import scipy
import scipy.ndimage
import scipy.stats
import matplotlib.pyplot as plt
from aux_functions import add_chr_to_bedtool, remove_chr_bedtool, add_chr_to_list
import pybedtools as pbt
import numpy as np

arr = np.asarray

def get_closest_atac_peaks_OLD(places, ALL_ATAC_PEAKS, k=10):
    z = places.sort().closest(ALL_ATAC_PEAKS, k=k, d=True)
    atac_peakset = []
    for i in z:
        atac_reg = tuple(i[5:8])
        atac_peakset.append(atac_reg)
    return pbt.BedTool(list(atac_peakset))

def get_closest_atac_peaks(places, ALL_ATAC_PEAKS, k=10):
    z = places.sort().closest(ALL_ATAC_PEAKS, k=k, d=True).to_dataframe(header=None)
    return z


def get_closest_atac_peak_dict(places, ALL_ATAC_PEAKS, k=10):
    z = places.sort().closest(ALL_ATAC_PEAKS, k=k, d=True)
    atac_peakset = []
    for i in z:
        atac_reg = tuple(i[5:8])
        atac_peakset.append(atac_reg)
    return pbt.BedTool(list(atac_peakset))


def filter_bad_genes(genelist):
    newgenelist = []
    for x in genelist:
        if ('Gm' in x) or ('Rik' in x) or ('-ps' in x) or ("-mt" in x) or ("Rp" in x):
            continue
        else:
            newgenelist.append(x)
    return newgenelist

from collections import defaultdict
def aggregate_binding_at_nearby_foxp3_peaks_OLD(chip_dict, deg_dict, foxp3_peaks_yuri, 
                            mm10_genes, PARSED_CHROMS, K=5, vmax=100, useControl=True, delta=10_000):
    if useControl:
        control_foxp3_peaks = add_chr_to_bedtool(foxp3_peaks_yuri)
        collapse_control_val_dict = {}
        control_val_dict = {}
        for name, bw in chip_dict.items():
            vs = fetch_from_bw(bw, control_foxp3_peaks, [], delta, bins=1_000)
            collapse_vs = np.nanmean(vs, axis=0)
            collapse_vs = scipy.ndimage.gaussian_filter1d(collapse_vs, sigma=40)
            collapse_control_val_dict[name] = collapse_vs
            control_val_dict[name] = vs

    n = len(deg_dict)
    v_dict = defaultdict(dict)
    fig, axs = initialize_subplots(n, 2, fgsz=4)
    for c, (cond, genelist) in enumerate(deg_dict.items()):
        ax = axs[c]
        if len(genelist)==0:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(cond)
            continue
        pos_places, neg_places = genes_to_start(mm10_genes, genelist, PARSED_CHROMS)
        places = remove_chr_bedtool(pbt.BedTool(pos_places + neg_places))
        atac_peaks = add_chr_to_bedtool(get_closest_atac_peaks(places, foxp3_peaks_yuri, k=K))

        pvals = []
        all_vs = []
        for name, bw in chip_dict.items():
            vs = fetch_from_bw(bw, atac_peaks, [], delta, bins=1_000)
            all_vs.append(vs)
            collapse_vs = np.nanmean(vs, axis=0)
            collapse_vs = scipy.ndimage.gaussian_filter1d(collapse_vs, sigma=40)
            ax.plot(collapse_vs, color='blue', alpha=.2)

            if useControl:
                control_collapse_vs = collapse_control_val_dict[name]
                ax.plot(control_collapse_vs, color='red', alpha=.2)
                center_control = control_val_dict[name][:, delta]
                center_deg = vs[:, delta]
                pval = scipy.stats.ranksums(center_deg, center_control)[1]
                pvals.append(pval.round(2))
            v_dict[cond][name] = vs
        ax.set_title(f'{cond}\n{pvals}')
        ax.set_ylim([0, vmax])
        blue_line = mlines.Line2D([], [], color='blue', marker='_', markersize=15, label='DEGs')
        if useControl:
            red_line = mlines.Line2D([], [], color='red', marker='_', markersize=15, label='Control')
            ax.legend(handles=[blue_line, red_line])
    for ax in axs[c+1:]:
        ax.set_xticks([])
        ax.set_yticks([])
    plt.tight_layout()
    return v_dict

def aggregate_binding_at_TSS_by_K_nearest_peaks(tss_df, foxp3_peaks_yuri, bw, K=5, delta=5_000):
    tss_bedtool = pbt.BedTool.from_dataframe(tss_df)
    atac_peaks = get_closest_atac_peaks(tss_bedtool, foxp3_peaks_yuri, k=K)
    atac_peak_bedtool = add_chr_to_bedtool(pbt.BedTool.from_dataframe(atac_peaks.iloc[:, 4:]))
    vs = fetch_from_bw(bw, atac_peak_bedtool, [], delta, bins=1_000)
    return vs,  atac_peaks

def aggregate_binding_at_peakset(peak_bedtool, bw, delta=5_000, bins=1000):
    vs = fetch_from_bw(bw, peak_bedtool, [], delta, bins=bins)
    return vs, peak_bedtool.to_dataframe(header=None)


def fetch_from_bw(bw, pos_places, neg_places, delta, **kwargs):
    vs = []
    if pos_places:
        chrom, s, e = unzip_for_bw(pos_places)[:3]
        s = (arr(s).astype(int)+arr(e).astype(int))//2
        vs_pos = bw.stackup(chrom, s-delta, s+delta, **kwargs)
        vs.append(vs_pos)
        
    if neg_places:
        chrom, s, e = unzip_for_bw(neg_places)[:3]
        s = (arr(s).astype(int)+arr(e).astype(int))//2
        vs_neg = bw.stackup(chrom, s-delta, s+delta, **kwargs)
        vs_neg = vs_neg[::-1]
        vs.append(vs_neg)
    vs = np.concatenate(vs, axis=0)
    return vs


from aux_functions import get_col
def get_average_at_bedtool(bedtool):
    chroms = get_col(bedtool, 0)
    starts = get_col(bedtool, 1).astype(float)
    ends = get_col(bedtool, 2).astype(float)
    return pbt.BedTool(list(zip(chroms, starts, ends)))
    

def aggregate_binding_at_sites(chip_dict, place_dict):
    v_dict = defaultdict(dict)
    for c, (cond, places) in enumerate(place_dict.items()):
        places = pbt.BedTool(places)
        all_vs = []
        for name, bw in chip_dict.items():
            vs = fetch_from_bw(bw, places, [], 10_000, bins=1000)
            all_vs.append(vs)
            collapse_vs = np.nanmean(vs, axis=0)
            collapse_vs = scipy.ndimage.gaussian_filter1d(collapse_vs, sigma=40)
            v_dict[cond][name] = vs
    return v_dict


def aggregate_chip_values_by_genes(values, values_metadata, index = [500], genesoi=None):
    # smoothed_values = scipy.ndimage.gaussian_filter1d(values, sigma=10)
    smoothed_values = values
    rows = []
    if genesoi is None:
        genesoi = values_metadata['gene'].unique()
    gene_values = []
    all_values = []
    for u in genesoi:
        indsoi = values_metadata['gene'] == u
        gene_value = smoothed_values[indsoi][:, index].mean()
        rows.append([u, gene_value])
        gene_values.append(smoothed_values[indsoi].mean(axis=0))
        all_values.append(smoothed_values[indsoi])
    return pd.DataFrame(rows, columns=['gene', 'value']), gene_values, all_values



import math
def initialize_subplots(n, rows, fgsz=4):
    cols = math.ceil(n/rows)
    fig, axs = plt.subplots(rows, cols, figsize=(cols*fgsz, rows*fgsz))
    axs = np.ravel(axs)
    return fig, axs

def unzip_for_bw(places):
    return list(zip(*places))


def genes_to_start(mm10_gff, genelist, PARSED_CHROMS, delta = 0):
    pos_places = []
    neg_places = []
    resultdict = {}
    names = []
    for i in mm10_gff:
        genename = i[-1].split('Name=')[1].split(';')[0]
        if genename in genelist:
            chrom, s, e = i[0], i[3], i[4]
            chrom = add_chr_to_list([chrom])[0]
            
            if chrom not in PARSED_CHROMS:
                continue
            
            strand = i[6]
            if strand == '+':
                pos_places.append([chrom, int(s)-delta, int(s)+delta, genename, strand])
            elif strand == '-':
                neg_places.append([chrom, int(e)-delta, int(e)+delta, genename, strand])
    return pos_places, neg_places



import pandas as pd
def gene_to_closest_atac_peaks(tss_df, atac_bedtool):
    tss_bedtool = pbt.BedTool.from_dataframe(tss_df.reset_index()[['chrom', 'start', 'end', 'gene_name']])
    gene_to_closest_atac_peaks_by_K = {}
    for k in [5, 10, 15]:
        gene_to_closest_atac_peaks = defaultdict(list)

        z = tss_bedtool.sort().closest(atac_bedtool, k=k, d=True, t='first')
        for row in z:
            genename = row[3]
            atac_peak = '_'.join(row[4:7])
            gene_to_closest_atac_peaks[genename].append(atac_peak)
            assert len(gene_to_closest_atac_peaks[genename]) <= k
        gene_to_closest_atac_peaks_by_K[k] = gene_to_closest_atac_peaks

    gene_to_closest_atac_peak_df_dict = {}
    for k, gene_to_closest_atac_peaks in gene_to_closest_atac_peaks_by_K.items():
        cols = []
        rows = []
        for i, row in gene_to_closest_atac_peaks.items():
            if len(row) < 2:
                continue
            cols.append(i)
            rows.append(row)
        gene_to_closest_atac_peak_df = pd.DataFrame(rows).T
        gene_to_closest_atac_peak_df.columns = cols
        gene_to_closest_atac_peak_df_dict[k] = gene_to_closest_atac_peak_df    
    return gene_to_closest_atac_peaks_by_K, gene_to_closest_atac_peak_df_dict


import pandas as pd
def gene_to_nearby_atac_peaks(tss_df, atac_bedtool):
    tss_bedtool = pbt.BedTool.from_dataframe(tss_df.reset_index()[['chrom', 'start', 'end', 'gene_name']])
    gene_to_nearby_atac_peaks_by_K = {}
    for k in [10_000, 40_000, 100_000]:
        gene_to_closest_atac_peaks = defaultdict(list)
        z = tss_bedtool.slop(b=k, g='./annotations/chromsizes').intersect(atac_bedtool, wo=True)
        for row in z:
            genename = row[3]
            atac_peak = '_'.join(row[4:7])
            gene_to_closest_atac_peaks[genename].append(atac_peak)
            assert len(gene_to_closest_atac_peaks[genename]) <= k
        gene_to_nearby_atac_peaks_by_K[k] = gene_to_closest_atac_peaks
    return gene_to_nearby_atac_peaks_by_K

def get_atac_peaks(deg_dict, gene_to_peak_dict, cond):
    geneset = deg_dict[cond]
    genenames = [x for x in geneset if x in gene_to_peak_dict.keys()]
    atac_peaks_of_interest = []
    for x in genenames:
        atac_peaks_of_interest.extend(gene_to_peak_dict[x] )
    return np.ravel(atac_peaks_of_interest)



import seaborn as sns
def chip_enrichment_by_geneset(deg_dict, gene_to_peak_dict, atac_metadata_df, cond1, cond2, xlabel, ylabel, title="", ax = None):
    enrichment_df = pd.DataFrame()
    for label, (geneset) in deg_dict.items():
        genenames = [x for x in geneset if x in gene_to_peak_dict.keys()]
        atac_peaks_of_interest = []
        for x in genenames:
            atac_peaks_of_interest.extend(gene_to_peak_dict[x] )
        deg_subdf = atac_metadata_df.loc[np.ravel(atac_peaks_of_interest)]
        n_deg_peaks = len(deg_subdf)
        n_total_peaks = len(atac_metadata_df)
        p_deg_peak = n_deg_peaks/n_total_peaks
        baseline = atac_metadata_df.mean(axis=0)
        p_both = deg_subdf.sum(axis=0)/n_total_peaks

        enrichments = p_both/(baseline*p_deg_peak)
        enrichment_df[label] = enrichments
        enrichment_df.index = atac_metadata_df.columns


    notnan_df = enrichment_df.loc[:, ~enrichment_df.isna().any(axis=0)]
    log_enrich_df = np.log2(notnan_df)
    log_enrich_df[np.isinf(log_enrich_df)]=-10

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(4, 4))
    x = log_enrich_df[cond1]
    y = log_enrich_df[cond2]
    bad = x.index.str.contains("_vs_") | (x < -4) | (y < -4)
    x = x[~bad]
    y = y[~bad]
    z = np.linspace(-3, 3)
    plt.sca(ax)
    plt.scatter(x[~bad], y[~bad], s=8)
    colors = sns.color_palette('tab20')[2::2]
    for c, index in enumerate(['Satb1', 'Ep300', 'Foxp3', 'Foxp1', 
                            'H3K27me3', 'Med1', 'Relap', 'IRF4', 'TCF1']):
        indsoi = x.index.str.contains(index)
        plt.scatter(x[indsoi], y[indsoi], s=18, color=colors[c], label=index)
    plt.plot(z, z, 'k', zorder=-1)
    plt.xlim([-5, 5])
    plt.ylim([-5, 5])
    plt.xlabel(f"{xlabel}")
    plt.ylabel(f"{ylabel}")
    plt.title(f"{title}")
    plt.tight_layout()    
    plt.suptitle(f"{title}", va='bottom')
    return log_enrich_df
