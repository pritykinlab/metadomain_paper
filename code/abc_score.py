from collections import defaultdict
from aux_functions import samesize_nonan_test, remove_chr_bedtool
import pybedtools as pbt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotting_functions import init_subplots_exact, mm

arr = np.asarray


def shuffle_diagonals(matrix, method='intra'):
    n = len(matrix)
    # Create a copy to avoid modifying the original matrix
    shuffled_matrix = np.copy(matrix)
    
    # Handle diagonals below and including the main diagonal
    if method == 'intra':
        for offset in range(0, n):
            # Extract diagonal elements
            diag_indices = np.diag_indices(n - offset)
            diagonal = shuffled_matrix[diag_indices[0], diag_indices[1] + offset]
            
            # Shuffle the diagonal elements
            diagonal = np.random.permutation(diagonal)
            
            # Assign the shuffled elements back to the matrix
            shuffled_matrix[diag_indices[0], diag_indices[1] + offset] = diagonal
        shuffled_matrix = np.triu(shuffled_matrix) + np.triu(shuffled_matrix, 1).T
    elif method == 'inter':
        for offset in range(0, n):
            # Extract diagonal elements
            shuffled_matrix[offset, :] = np.random.permutation(shuffled_matrix[offset, :])
    else:
        raise ValueError("Invalid randomization method. Choose either 'intra' or 'inter'.")
    return shuffled_matrix

def compute_ABC_scores(coolfile, h3k7ac_vec, chrom_to_start, chrom_to_end, parsed_chroms, use_chr=False, res=50_000, randomize_hic=None):
    intra_activity_dfs = defaultdict(list)
    inter_activity_dfs = defaultdict(list)
    for chrom1 in parsed_chroms:
        for chrom2 in parsed_chroms:
            if use_chr:
                contact = coolfile.matrix().fetch('chr' + chrom1, 'chr' + chrom2)
            else:
                contact = coolfile.matrix().fetch(chrom1, chrom2)
            if randomize_hic == 'fullrandom':
                contact = np.random.permutation(contact)
            elif randomize_hic == 'diagrandom':
                if chrom1 == chrom2:
                    contact = shuffle_diagonals(contact, method='intra')
                else:
                    contact = shuffle_diagonals(contact, method='inter')
            activity = np.ravel(h3k7ac_vec)[chrom_to_start[chrom2]:chrom_to_end[chrom2]]
            if chrom1 == chrom2:
                activity_df = pd.DataFrame()
                for cutoff in [10_000, 50_000, 100_000, 200_000, 500_000, 1e6, 2e6, 5e6, 10e6, 20e6, 40e6, 'full_intra']:
                    z = contact.copy()
                    z[np.isnan(z)] = 0
                    if cutoff != 'full_intra':
                        k = cutoff//res
                        z = np.triu(z, k=-k+1)
                        z = np.tril(z, k=k-1)
                    activities = z@(activity)
                    activity_df[f'{cutoff}'] = activities
                intra_activity_dfs[chrom1].append(activity_df)
            else:
                activity_df = pd.DataFrame()
                z = contact.copy()
                z[np.isnan(z)] = 0
                activities = z@(activity)
                activity_df[f'{chrom2}'] = activities
                inter_activity_dfs[chrom1].append(activity_df)
        print("Done with", chrom1)
    return intra_activity_dfs, inter_activity_dfs

def plot_abc(abc_dfs, tss_5kb_ind_df, peak_gene_dataframe, geneLengths, parsed_chroms, ylims = [0, .8]):
    rows = []
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 20*mm), dpi=200)
    for chrom in parsed_chroms:
        if chrom == 'X':
            continue
        genesoi = tss_5kb_ind_df[tss_5kb_ind_df['chr'] == chrom].index.intersection(
            peak_gene_dataframe['name'][peak_gene_dataframe['itemRgb'] > 4].values    
        ).intersection(geneLengths.index)
        idx = (tss_5kb_ind_df.loc[genesoi]['start'].astype(int).values/5000).astype(int)
        yvals = peak_gene_dataframe.set_index('name').loc[genesoi]['itemRgb']
        yvals = (1+yvals[~yvals.index.duplicated(keep='first')])
        yvals = yvals/geneLengths.loc[yvals.index]
        yvals = np.log2(yvals)

        rs = [] 
        abc_df = (abc_dfs[chrom][0].T - abc_dfs[chrom][0]['10000']).T
        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra']
        for col in cols:
            r = samesize_nonan_test(abc_df.iloc[idx][col].values, yvals)[0]
            rs.append(r)    
        cols[-1] = 'All'
        
        rows.append(rs)
        cols[:-1] = (arr(cols[:-1]).astype(float)/1e3).astype(int)
        plt.plot(cols, rs, marker='o', label='chr' + chrom, markersize=2)
        plt.xticks(rotation=30);
    plt.legend(bbox_to_anchor=(1, 1), ncol=2, fontsize=4, frameon=False)
    plt.ylim(ylims)
    plt.xlabel('Distance for ABC')
    plt.ylabel("Correlation")
    plt.title("Gene expression-ABC correlation, by distance")
    fig.savefig('./plots/ABC/ABC_gene_correlation_cumulative.pdf', bbox_inches = 'tight')    
    vals = pd.DataFrame(rows, columns = cols)
    return vals

import seaborn as sns
from aux_functions import nonan_test
from plot_pvals import add_stat_annotation_boxplot_no_hue, add_stat_annotation_boxplot_with_hue
def plot_abc2(abc_dfs, tss_5kb_ind_df, peak_gene_dataframe, geneLengths, parsed_chroms, ylims = [0, .8]):
    rows = []
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 20*mm), dpi=200)
    for chrom in parsed_chroms:
        if chrom == 'X':
            continue
        genesoi = tss_5kb_ind_df[tss_5kb_ind_df['chr'] == chrom].index.intersection(
            peak_gene_dataframe['name'][peak_gene_dataframe['itemRgb'] > 4].values    
        ).intersection(geneLengths.index)
        idx = (tss_5kb_ind_df.loc[genesoi]['start'].astype(int).values/5000).astype(int)
        yvals = peak_gene_dataframe.set_index('name').loc[genesoi]['itemRgb']
        yvals = (1+yvals[~yvals.index.duplicated(keep='first')])
        yvals = yvals/geneLengths.loc[yvals.index]
        yvals = np.log2(yvals)

        rs = [] 
        abc_df = (abc_dfs[chrom][0].T - abc_dfs[chrom][0]['10000']).T
        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra']
        col2name = {
            '50000': '50kb',
            '100000': '100kb',
            '1000000.0': '1Mb',
            '10000000.0': '10Mb',
            'full_intra': 'Whole\nChrom'
        }
        for col in cols:
            r = samesize_nonan_test(abc_df.iloc[idx][col].values, yvals)[0] 
            rows.append([r, chrom, col2name[col]])

    data = pd.DataFrame(rows, columns = ['r', 'chrom', 'distance'])
    sns.boxplot(data = data, x = 'distance', y = 'r', showfliers=False)
    sns.stripplot(data = data, x = 'distance', y = 'r', s=2, color='black')
    stat_pairs = []
    for i, col in enumerate(cols[:-1]):
        stat_pairs.append([col2name[col], col2name[cols[i+1]]])
    add_stat_annotation_boxplot_no_hue(plt.gca(), data, 'distance', 'r', [col2name[x] for x in cols], 
                                       stat_pairs, ymax=.37, delta = .18, h = 0.01,
                                       log=False)
    plt.ylabel("Correlation")
    plt.ylim(ylims)
    plt.xlabel('Distance for ABC')
    plt.ylabel("Correlation")
    plt.title("Gene expression-ABC correlation, by distance")
    fig.savefig('./plots/ABC/ABC_gene_correlation_cumulative.pdf', bbox_inches = 'tight')    
    # vals = pd.DataFrame(rows, columns = cols)
    # return vals


def plot_abc2(abc_dfs, tss_5kb_ind_df, peak_gene_dataframe, geneLengths, parsed_chroms, ylims = [0, .8]):
    rows = []
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 20*mm), dpi=200)
    for chrom in parsed_chroms:
        if chrom == 'X':
            continue
        genesoi = tss_5kb_ind_df[tss_5kb_ind_df['chr'] == chrom].index.intersection(
            peak_gene_dataframe['name'][peak_gene_dataframe['itemRgb'] > 4].values    
        ).intersection(geneLengths.index)
        idx = (tss_5kb_ind_df.loc[genesoi]['start'].astype(int).values/5000).astype(int)
        yvals = peak_gene_dataframe.set_index('name').loc[genesoi]['itemRgb']
        yvals = (1+yvals[~yvals.index.duplicated(keep='first')])
        yvals = yvals/geneLengths.loc[yvals.index]
        yvals = np.log2(yvals)

        rs = [] 
        abc_df = (abc_dfs[chrom][0].T - abc_dfs[chrom][0]['10000']).T
        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra']
        col2name = {
            '50000': '50kb',
            '100000': '100kb',
            '1000000.0': '1Mb',
            '10000000.0': '10Mb',
            'full_intra': 'Whole\nChrom'
        }
        for col in cols:
            r = samesize_nonan_test(abc_df.iloc[idx][col].values, yvals)[0] 
            rows.append([r, chrom, col2name[col]])

    data = pd.DataFrame(rows, columns = ['r', 'chrom', 'distance'])
    sns.boxplot(data = data, x = 'distance', y = 'r', showfliers=False)
    sns.stripplot(data = data, x = 'distance', y = 'r', s=2, color='black')
    stat_pairs = []
    for i, col in enumerate(cols[:-1]):
        stat_pairs.append([col2name[col], col2name[cols[i+1]]])
    add_stat_annotation_boxplot_no_hue(plt.gca(), data, 'distance', 'r', [col2name[x] for x in cols], 
                                       stat_pairs, ymax=.37, delta = .18, h = 0.01,
                                       log=False)
    plt.ylabel("Correlation")
    plt.ylim(ylims)
    plt.xlabel('Distance for ABC')
    plt.ylabel("Correlation")
    plt.title("Gene expression-ABC correlation, by distance")
    fig.savefig('./plots/ABC/ABC_gene_correlation_cumulative.pdf', bbox_inches = 'tight')    
    # vals = pd.DataFrame(rows, columns = cols)
    # return vals


def plot_abc3(abc_dfs_real, abc_dfs_control, tss_5kb_ind_df, peak_gene_dataframe, geneLengths, parsed_chroms, ylims = [0, .8]):
    rows = []
    for chrom in parsed_chroms:
        if chrom == 'X':
            continue
        genesoi = tss_5kb_ind_df[tss_5kb_ind_df['chr'] == chrom].index.intersection(
            peak_gene_dataframe['name'][peak_gene_dataframe['itemRgb'] > 4].values    
        ).intersection(geneLengths.index)
        idx = (tss_5kb_ind_df.loc[genesoi]['start'].astype(int).values/5000).astype(int)
        yvals = peak_gene_dataframe.set_index('name').loc[genesoi]['itemRgb']
        yvals = (1+yvals[~yvals.index.duplicated(keep='first')])
        yvals = yvals/geneLengths.loc[yvals.index]
        yvals = np.log2(yvals)

        rs = [] 
        abc_df = (abc_dfs_real[chrom][0].T - abc_dfs_real[chrom][0]['10000']).T
        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra']
        col2name = {
            '50000': '50kb',
            '100000': '100kb',
            '1000000.0': '1Mb',
            '10000000.0': '10Mb',
            'full_intra': 'Whole\nChrom'
        }
        for col in cols:
            r = samesize_nonan_test(abc_df.iloc[idx][col].values, yvals)[0] 
            rows.append([r, chrom, col2name[col], 'real'])

        abc_df = (abc_dfs_control[chrom][0].T - abc_dfs_control[chrom][0]['10000']).T
        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra']
        col2name = {
            '50000': '50kb',
            '100000': '100kb',
            '1000000.0': '1Mb',
            '10000000.0': '10Mb',
            'full_intra': 'Whole\nChrom'
        }
        for col in cols:
            r = samesize_nonan_test(abc_df.iloc[idx][col].values, yvals)[0] 
            rows.append([r, chrom, col2name[col], 'control'])

        

    data = pd.DataFrame(rows, columns = ['r', 'chrom', 'distance', 'control'])
    fig, axs = init_subplots_exact(1, 1, fgsz=(70*mm, 35*mm), dpi=200)
    sns.boxplot(data = data, x = 'distance', y = 'r', hue='control', showfliers=False,
                palette=['tab:blue', 'lightgray'])
    sns.stripplot(data = data, x = 'distance', y = 'r', hue='control', s=1.5, color='black',
                  dodge=True)
    stat_pairs = []
    for i, col in enumerate(cols[:-1]):
        if i == 0:
            stat_pairs.append([[col2name[col], 'real'], [col2name[cols[i]], 'control']])
        stat_pairs.append([[col2name[col], 'real'], [col2name[cols[i+1]], 'real']])
    add_stat_annotation_boxplot_with_hue(plt.gca(), data, 'distance', 'r', 'control', 
                                         [col2name[x] for x in cols], 
                                         ['real', 'control'],
                                         stat_pairs, ymax=.34, delta = .18, 
                                         yoff_method='None')
    plt.ylabel("Correlation")
    plt.ylim(ylims)
    plt.xlabel('Distance for ABC')
    plt.ylabel("Correlation")
    plt.title("Gene expression-ABC correlation, by distance")
    plt.legend(bbox_to_anchor=(1, 1))
    return fig
    # vals = pd.DataFrame(rows, columns = cols)
    # return vals




def plot_abc4(abc_dfs_real, abc_dfs_control, inter_abc_dfs_real, inter_abc_dfs_control, 
              all_genesoi, tss_5kb_ind_df, peak_gene_dataframe, geneLengths, parsed_chroms, ylims = [0, .8]):
    rows = []
    all_genesoi = np.asarray(all_genesoi)
    for chrom in parsed_chroms:
        if chrom == 'X':
            continue

        genesoi = all_genesoi[tss_5kb_ind_df.loc[all_genesoi, 'chr'] ==  chrom]
        idx = (tss_5kb_ind_df.loc[genesoi]['start'].astype(int).values/5000).astype(int)
        yvals = peak_gene_dataframe.set_index('name').loc[genesoi]['itemRgb']
        yvals = (1+yvals[~yvals.index.duplicated(keep='first')])
        yvals = yvals/geneLengths.loc[yvals.index]
        yvals = np.log2(yvals)

        rs = [] 

        abc_df = abc_dfs_real[chrom][0]

        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra', 'inter']
        col2name = {
            '50000': '50kb',
            '100000': '100kb',
            '1000000.0': '1Mb',
            '10000000.0': '10Mb',
            'full_intra': 'Whole\nChrom'
        }
        for col in cols:
            if col != 'inter':
                abc_vals = abc_df.iloc[idx][col] - abc_df.iloc[idx]['10000']
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, col2name[col], 'real'])
            else:
                abc_vals = abc_df.iloc[idx]['full_intra'].values - abc_df.iloc[idx]['10000']
                abc_vals = abc_vals + pd.concat(inter_abc_dfs_real[chrom], axis=1).sum(axis=1).iloc[idx]
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, 'Whole\nGenome', 'real'])
                

        abc_df = (abc_dfs_control[chrom][0])
        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra', 'inter']
        col2name = {
            '50000': '50kb',
            '100000': '100kb',
            '1000000.0': '1Mb',
            '10000000.0': '10Mb',
            'full_intra': 'Whole\nChrom',
            'inter' : 'Whole\nGenome',
        }
        for col in cols:
            if col != 'inter':
                abc_vals = abc_df.iloc[idx][col] - abc_df.iloc[idx]['10000']
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, col2name[col], 'control'])
            else:
                abc_vals = abc_df.iloc[idx]['full_intra'].values - abc_df.iloc[idx]['10000']
                abc_vals = abc_vals + pd.concat(inter_abc_dfs_control[chrom], axis=1).sum(axis=1).iloc[idx]
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, 'Whole\nGenome', 'control'])


    data = pd.DataFrame(rows, columns = ['r', 'chrom', 'distance', 'control'])
    fig, axs = init_subplots_exact(1, 1, fgsz=(70*mm, 35*mm), dpi=200)
    sns.boxplot(data = data, x = 'distance', y = 'r', hue='control', showfliers=False,
                palette=['tab:blue', 'lightgray'])
    sns.stripplot(data = data, x = 'distance', y = 'r', hue='control', s=1.5, color='black',
                  dodge=True)
    stat_pairs = []
    for i, col in enumerate(cols[:-1]):
        if i == 0:
            stat_pairs.append([[col2name[col], 'real'], [col2name[cols[i]], 'control']])
        stat_pairs.append([[col2name[col], 'real'], [col2name[cols[i+1]], 'real']])
    add_stat_annotation_boxplot_with_hue(plt.gca(), data, 'distance', 'r', 'control', 
                                         [col2name[x] for x in cols], 
                                         ['real', 'control'],
                                         stat_pairs, ymax=.34, delta = .18, 
                                         yoff_method='None')
    plt.ylabel("Correlation")
    plt.ylim(ylims)
    plt.yticks([0, .2, .4, .6, .8])
    plt.xlabel('Distance for ABC')
    plt.ylabel("Correlation")
    plt.title("RPKM-cABC correlation")
    plt.legend(bbox_to_anchor=(1, 1))
    return fig
    # vals = pd.DataFrame(rows, columns = cols)
    # return vals



def plot_abc_lfc(abc_dfs_real_treg, abc_dfs_control_treg, inter_abc_dfs_real_treg, inter_abc_dfs_control_treg, 
                 abc_dfs_real_tcon, abc_dfs_control_tcon, inter_abc_dfs_real_tcon, inter_abc_dfs_control_tcon, 
              all_genesoi, tss_5kb_ind_df, peak_gene_dataframe, parsed_chroms, ylims = [0, .8]):
    rows = []
    all_genesoi = np.asarray(all_genesoi)
    for chrom in parsed_chroms:
        if chrom == 'X':
            continue

        genesoi = all_genesoi[tss_5kb_ind_df.loc[all_genesoi, 'chr'] ==  chrom]
        idx = (tss_5kb_ind_df.loc[genesoi]['start'].astype(int).values/5000).astype(int)
        yvals = peak_gene_dataframe.set_index('name').loc[genesoi]['thickStart']

        rs = [] 

        abc_df_treg = abc_dfs_real_treg[chrom][0]
        abc_df_tcon = abc_dfs_real_tcon[chrom][0]

        col2name = {
            '50000': '50kb',
            '100000': '100kb',
            '1000000.0': '1Mb',
            '10000000.0': '10Mb',
            'full_intra': 'Whole\nChrom',
            'inter': 'Whole\nGenome',
        }
        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra', 'inter']
        for col in cols:
            if col != 'inter':
                abc_vals_treg = abc_df_treg.iloc[idx][col] - abc_df_treg.iloc[idx]['10000']
                abc_vals_tcon = abc_df_tcon.iloc[idx][col] - abc_df_tcon.iloc[idx]['10000']
                abc_vals = np.log2(abc_vals_treg / abc_vals_tcon)
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, col2name[col], 'real'])
            else:
                abc_vals_treg = abc_df_treg.iloc[idx]['full_intra'] - abc_df_treg.iloc[idx]['10000'] + pd.concat(inter_abc_dfs_real_treg[chrom], axis=1).sum(axis=1).iloc[idx]
                abc_vals_tcon = abc_df_tcon.iloc[idx]['full_intra'] - abc_df_tcon.iloc[idx]['10000'] + pd.concat(inter_abc_dfs_real_tcon[chrom], axis=1).sum(axis=1).iloc[idx]
                abc_vals = np.log2(abc_vals_treg / abc_vals_tcon)
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, 'Whole\nGenome', 'real'])
                

        abc_df_treg = abc_dfs_control_treg[chrom][0]
        abc_df_tcon = abc_dfs_control_tcon[chrom][0]

        cols = ['50000', '100000', '1000000.0', '10000000.0', 'full_intra', 'inter']
        for col in cols:
            if col != 'inter':
                abc_vals_treg = abc_df_treg.iloc[idx][col] - abc_df_treg.iloc[idx]['10000']
                abc_vals_tcon = abc_df_tcon.iloc[idx][col] - abc_df_tcon.iloc[idx]['10000']
                abc_vals = np.log2(abc_vals_treg / abc_vals_tcon)
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, col2name[col], 'control'])
            else:
                abc_vals_treg = abc_df_treg.iloc[idx]['full_intra'] - abc_df_treg.iloc[idx]['10000'] + pd.concat(inter_abc_dfs_control_treg[chrom], axis=1).sum(axis=1).iloc[idx]
                abc_vals_tcon = abc_df_tcon.iloc[idx]['full_intra'] - abc_df_tcon.iloc[idx]['10000'] + pd.concat(inter_abc_dfs_control_tcon[chrom], axis=1).sum(axis=1).iloc[idx]
                abc_vals = np.log2(abc_vals_treg / abc_vals_tcon)
                r = samesize_nonan_test(abc_vals.values, yvals)[0] 
                rows.append([r, chrom, 'Whole\nGenome', 'control'])


    data = pd.DataFrame(rows, columns = ['r', 'chrom', 'distance', 'control'])
    fig, axs = init_subplots_exact(1, 1, fgsz=(70*mm, 35*mm), dpi=200)
    sns.boxplot(data = data, x = 'distance', y = 'r', hue='control', showfliers=False,
                palette=['tab:blue', 'lightgray'])
    sns.stripplot(data = data, x = 'distance', y = 'r', hue='control', s=1.5, color='black',
                  dodge=True)
    stat_pairs = []
    for i, col in enumerate(cols[:-2]):
        if i == 0:
            stat_pairs.append([[col2name[col], 'real'], [col2name[cols[i]], 'control']])
        stat_pairs.append([[col2name[col], 'real'], [col2name[cols[i+1]], 'real']])
        stat_pairs.append([[col2name[col], 'real'], [col2name[cols[i+2]], 'real']])
    add_stat_annotation_boxplot_with_hue(plt.gca(), data, 'distance', 'r', 'control', 
                                         [col2name[x] for x in cols], 
                                         ['real', 'control'],
                                         stat_pairs, ymax=.45, delta = .11, 
                                         yoff_method='None')
    plt.ylabel("Correlation")
    plt.ylim(ylims)
    plt.yticks([0, .2, .4, .6, .8, 1])
    plt.xlabel('Distance for ABC')
    plt.ylabel("Correlation")
    plt.title("RPKM-cABC correlation")
    plt.legend(bbox_to_anchor=(1, 1))
    return fig



def make_tss_ind_df(all_ind_to_region_5kb, my_tss_df):
    rows = []
    seen = set()
    for i in pbt.BedTool(all_ind_to_region_5kb).intersect(remove_chr_bedtool(pbt.BedTool.from_dataframe(my_tss_df)), wo=True):
        name = i[6]
        region = i[:3]
        if name in seen:
            continue
        else:
            rows.append(region + [name])        
            seen.add(name)
    tss_5kb_ind_df = pd.DataFrame(rows, columns = ['chr', 'start', 'end', 'gene_name']).set_index('gene_name')    
    return tss_5kb_ind_df


from aux_functions import get_col
def make_ternary_df(intra_activity_dfs, inter_activity_dfs, all_ind_to_region, tss_df, TSS_filter=True, genesets=[]):
    vs_local = []
    vs_far = []
    vs_inter = []
    for chrom in inter_activity_dfs:
        local_intra = pd.concat(intra_activity_dfs[chrom], axis=1)['2000000.0'] - pd.concat(intra_activity_dfs[chrom], axis=1)['10000']
        full_intra = pd.concat(intra_activity_dfs[chrom], axis=1)['full_intra'] - pd.concat(intra_activity_dfs[chrom], axis=1)['10000']
        full_inter = pd.concat(inter_activity_dfs[chrom], axis=1).sum(axis=1)
        
        vs_local.extend(list((local_intra)/(full_intra + full_inter)))
        vs_far.extend(list((full_intra-local_intra)/(full_intra + full_inter)))
        vs_inter.extend(list(full_inter/(full_intra+full_inter)))

    
    vs_local_sum = []
    vs_far_sum = []
    vs_inter_sum = []
    for chrom in inter_activity_dfs:
        local_intra = pd.concat(intra_activity_dfs[chrom], axis=1)['2000000.0'] - pd.concat(intra_activity_dfs[chrom], axis=1)['10000']
        full_intra = pd.concat(intra_activity_dfs[chrom], axis=1)['full_intra'] - pd.concat(intra_activity_dfs[chrom], axis=1)['10000']
        full_inter = pd.concat(inter_activity_dfs[chrom], axis=1).sum(axis=1)
        
        vs_local_sum.extend(list(local_intra))
        vs_far_sum.extend(list(full_intra-local_intra))
        vs_inter_sum.extend(list(full_inter))
    ternary_df = pd.DataFrame([vs_local, vs_far, vs_inter, vs_local_sum, vs_far_sum, vs_inter_sum], 
                              index = ['Local', 'Far', 'Inter', 'Local_Sum', 'Far_Sum', 'Inter_Sum']).T

    
    ternary_df['chrom'] = [all_ind_to_region[x][0] for x in ternary_df.index]

    ternary_dfs = []
    for geneset in genesets:
        subdf = tss_df.loc[[x for x in geneset if x in tss_df.index]]
        tss_cols = get_col(pbt.BedTool(all_ind_to_region).intersect(pbt.BedTool.from_dataframe(subdf), c=True), -1).astype(int)

        sub_ternary_df = ternary_df.iloc[tss_cols>0].dropna()
        ternary_dfs.append(sub_ternary_df)
    return ternary_dfs