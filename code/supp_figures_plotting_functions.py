
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2

### Table S1
def make_table_s1():
    s = 0
    rows = []
    for rep in ['Tn.1', 'Tn.2', 'Tn.3', 'Treg.1', 'Treg.2', 'Treg.3']:
        total = pd.read_csv(f'/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/logs/OGYuri.{rep}.base_pairs_stats',
                    sep='\t', header=None).set_index(0).loc['total'].values[0]
        total_mapped = pd.read_csv(f'/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/logs/OGYuri.{rep}.base_pairs_stats',
                    sep='\t', header=None).set_index(0).loc['total_mapped'].values[0]
        
        total_nodups = pd.read_csv(f'/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/logs/OGYuri.{rep}.deduped_pair_stats',
                    sep='\t', header=None).set_index(0).loc['total_mapped'].values[0]
        
        cis = pd.read_csv(f'/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/logs/OGYuri.{rep}.deduped_pair_stats',
                    sep='\t', header=None).set_index(0).loc['cis'].values[0]
        
        trans = pd.read_csv(f'/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/logs/OGYuri.{rep}.deduped_pair_stats',
                    sep='\t', header=None).set_index(0).loc['trans'].values[0]
        s += total_nodups
        row = pd.Series([total, total_mapped, total_nodups, cis, trans],
                    index = ['total', 'total_properly_mapped', 'total_nodups', 'intra', 'inter'],
                    name = rep.replace("Tn", "Tcon"))
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv('./plots/paper/tables/table_s1.tsv', sep='\t')
    return table


### Fig S1
def make_hic_jesseye_readcount_barplot():
    data = []
    for file in ['/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_merged_cools/TReg.mcool',
                '/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_merged_cools/TCon.mcool']:
        name = file.split("/")[-1].split('.mcool')[0]
        sum = cooler.Cooler(file + "::/resolutions/250000").info['sum']
        data.append([file, name.lower().capitalize(), 'Jesse/Ye', sum, ])

    for file in ['/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_merged_cools/Rename_Tconv_all_no_chrM.mcool',
    '/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_merged_cools/Rename_Treg_all_no_chrM.mcool']:
            name = file.split("/")[-1].split('.mcool')[0]
            sum = cooler.Cooler(file + "::/resolutions/250000").info['sum']
            data.append([file, name.split("_")[1], 'Ours', sum, ])

    data = pd.DataFrame(data, columns = ['path', 'name', 'source', 'sum'])
    data['sum'] = data['sum']/1e6

    data = data.groupby(['name', 'source']).sum().sort_values('sum')    
    fig, axs = init_subplots_exact(1, 1, fgsz=(60*mm, 60*mm), dpi=100)
    data['sum'].plot.bar(ax=axs)
    plt.xticks(rotation=10)
    plt.xlabel("Dataset")
    plt.ylabel("M Reads")
    plt.title("M reads in dataset")
    fig = plt.gcf()
    return fig    

def make_venn_diagrams(all_loops_df):
    figs = []
    sets = []
    colsoi = all_loops_df.columns[all_loops_df.columns.str.contains("Treg_rep1") | all_loops_df.columns.str.contains("Treg_all")]
    for i in colsoi:
        sets.append(set(np.where(all_loops_df[i])[0]))
    
    colnames = [x.split("supplement_")[1] for x in colsoi]
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)
    venn2(sets, colnames, ax=axs)
    
    sets = []
    colsoi = all_loops_df.columns[all_loops_df.columns.str.contains("Tn_rep1") | all_loops_df.columns.str.contains("Tconv_all")]
    for i in colsoi:
        sets.append(set(np.where(all_loops_df[i])[0]))
    
    colnames = [x.split("supplement_")[1] for x in colsoi]
    figs.append(plt.gcf())
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)
    venn2(sets, colnames, ax=axs)
    
    sets = []
    colsoi = all_loops_df.columns[all_loops_df.columns.str.contains("Treg_all") 
                                  | all_loops_df.columns.str.contains("Tconv_all")
                                  | all_loops_df.columns.str.contains("Merged_all")]
    for i in colsoi:
        sets.append(set(np.where(all_loops_df[i])[0]))
    
    colnames = [x.split("supplement_")[1] for x in colsoi]
    figs.append(plt.gcf())
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)
    venn3(sets, colnames, ax=axs)
    figs.append(plt.gcf())
    return figs

import numpy as np
import seaborn as sns

def make_jaccard_clustermap(all_loops_df):
    intersection = all_loops_df.T @ all_loops_df
    
    v = all_loops_df.sum(axis=0)
    fraction_intersecting = intersection / (np.add.outer(v.values, v.values) - intersection)
    colsoi = fraction_intersecting.index.str.contains("rep")
    df = fraction_intersecting.loc[colsoi, colsoi]
    df.columns = df.columns.str.replace(".*T", 'T', regex=True)
    df.index = df.index.str.replace(".*T", 'T', regex=True)
    g = sns.clustermap(df, cmap='bwr', vmin=.1, vmax=.8, annot=True, zorder=3, figsize=(8, 8))
    return g.fig



import seaborn as sns
import matplotlib.pyplot as plt
from plotting_functions import mm

def make_fdr_histogram(init_dfs):
    init_df_oi = init_dfs['final_loop_call_for_paper_supplement_Merged_all'].copy()
    fdr_values = init_df_oi['FDR']
    
    fig, ax = plt.subplots(1, 1, figsize=(40*mm, 40*mm), dpi=200)
    sns.histplot(fdr_values, ax=ax)
    ax.set_title("Loop FDRs from ")
    ax.set_axisbelow(True)
    plt.axvspan(0, .02, color='black', alpha=.2)
    plt.axvline(.02, color='black', linestyle='--')
    # plt.xlim()
    return fig


import numpy as np
import random
from plotting_functions import *
from aux_functions import *

from process_loops_from_mustache import *
print("Imported plotting functions")
def plot_random_examples(init_dfs, cooldict_5kb):
    unmerged_loopdict = {
        'merged': loops_from_df(init_dfs['final_loop_call_for_paper_supplement_Merged_all']),
        'tcon': loops_from_df(init_dfs['final_loop_call_for_paper_supplement_Tconv_all']),
        'treg': loops_from_df(init_dfs['final_loop_call_for_paper_supplement_Treg_all']),
    }
    
    random.seed(0)
    
    init_df_oi = init_dfs['final_loop_call_for_paper_supplement_Merged_all'].copy()
    init_df_oi.index = loopset_to_granges(loops_from_df(init_df_oi))
    
    fdr_values = init_df_oi['FDR']
    digitized_fdr_values = np.digitize(fdr_values, [0, 1e-4, 1e-3, 1e-2, 5e-2])
    
    figs = []
    for u in np.unique(digitized_fdr_values):
        fdr_bool = digitized_fdr_values == u
        granges = random.choices(init_df_oi[fdr_bool].index, k=6)
        pmin, pmax = fdr_values[fdr_bool].min(), fdr_values[fdr_bool].max()
        pmin, pmax = map(get_pval_string_from_pval, [pmin, pmax])
        fig = plot_from_one_cooler_and_multiple_granges(
            cooldict_5kb, granges, 5000 * 15, loopdict=unmerged_loopdict, cond='merged',
            suptitle=f'P-value quantile {u}: {pmin} -- {pmax}',
            add_chr=True,
        )
        plt.gcf().savefig(f'./plots/paper/s2/quantile_{u}.pdf', bbox_inches='tight')
        # figs.append(fig)
    # return figs

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def make_manhattan_distance_plot(all_loops):
    loop_df = all_loops.to_dataframe(header=None).iloc[:, :6]
    loop_df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    mindists = []
    for chrom in loop_df['chrom1'].unique():
        subdf = loop_df[loop_df['chrom1'] == chrom]
        vx = np.add.outer(subdf['start1'].values, -subdf['start1'].values).astype(float)
        vy = np.add.outer(subdf['start1'].values, -subdf['start2'].values).astype(float)
        v = np.abs(vx) + np.abs(vy)
        v += np.diag(np.diag(v) * np.nan)
        mindist = np.nanmin(v, axis=1)
        mindists.extend(list(mindist))    
    manhattan_distance = np.array(mindists)
    fig, ax = plt.subplots(1, 1, figsize=(40*mm, 40*mm), dpi=200)
    sns.histplot(manhattan_distance, stat='count', log_scale=True, ax=ax, rasterized=True)
    ax.set_xlabel("Manhattan Distance (bp)")
    ax.set_title("Loop manhattan distance")
    ax.set_axisbelow(True)
    fig.savefig('./plots/paper/s2/loop_manhattan.pdf', bbox_inches='tight', dpi = 500)
    plt.show()


def deseq_ma_plot(all_loops):
    c = (np.sign(all_loops[7]) * (all_loops[11] < .05))
    fig, axs = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi = 150)
    plt.scatter(all_loops[6], all_loops[7].clip(-4, 4), s = .4,
            c = c, cmap='coolwarm')
    plt.xscale('log')
    plt.ylabel("Loop LFC")
    plt.xlabel("Loop Read Count")
    plt.text(1, 0, f'NS={(c==0).sum():,}', va='bottom', ha='right', transform=axs.transAxes,
            fontsize=6)
    plt.text(1, .08, f'Treg={(c==1).sum():,}', va='bottom', ha='right', transform=axs.transAxes, color='red',
            fontsize=6)
    plt.text(1, .16, f'Tcon={(c==-1).sum():,}', va='bottom', ha='right', transform=axs.transAxes, color='blue',
            fontsize=6)    
    return fig


import glob
import os
import pandas as pd
import pybedtools as pbt

def adjust_index_to_be_5kb(loop_df, delta=15_000):
    adj_index = []
    for i in loop_df.index:
        l = loop_grange_to_tuple(i)
        l1, l2 = l[:3], l[3:6]
        l1, l2 = map(make_int, [l1, l2])
        l1, l2 = map(lambda x: extend_l(x, -delta), [l1, l2])
        l1, l2 = map(lambda x: tuple_to_grange(*x), [l1, l2])        
        loop_grange = granges_to_loop_grange(l1, l2)
        adj_index.append(loop_grange)
    loop_df.index = adj_index
    return loop_df

def get_ancs_from_loops(loop_df):
    rows = []
    for x in loop_df.index:
        grange1, grange2 = x.split("|")
        anc1 = grange_to_tuple(grange1)
        anc2 = grange_to_tuple(grange2)
        rows.append(list(anc1) + list(anc2))
    anc_df = pd.DataFrame(rows, columns=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
    anc_df.index = loop_df.index
    return anc_df

import pandas as pd
import numpy as np
import networkx as nx

def compare_jesseye_loops_to_our_loops(loopdict):

    jesse_ye_loopdata = pd.read_csv('./jesse_ye_diff_loops/jesse_ye_diff_loops.csv')

    def liftover_df(df, lo):
        def liftover_coord(chromosome, start, end):
            lifted_start = lo.convert_coordinate(chromosome, start)
            lifted_end = lo.convert_coordinate(chromosome, end)
            if lifted_start and lifted_end:
                return lifted_start[0][1], lifted_end[0][1]
            else:
                return None, None

        # Apply liftover to each row
        new_coords_bin1 = df.apply(lambda row: liftover_coord(row['chromosome'], row['bin1_start'], row['bin1_end']), axis=1)
        df['bin1_start'], df['bin1_end'] = zip(*new_coords_bin1)

        new_coords_bin2 = df.apply(lambda row: liftover_coord(row['chromosome'], row['bin2_start'], row['bin2_end']), axis=1)
        df['bin2_start'], df['bin2_end'] = zip(*new_coords_bin2)
        return df

    def swap_bins(df):
        # Swap bin1 and bin2 if bin1_start is greater than bin2_start
        def swap_if_needed(row):
            if row['bin1_start'] > row['bin2_start']:
                row['bin1_start'], row['bin2_start'] = row['bin2_start'], row['bin1_start']
                row['bin1_end'], row['bin2_end'] = row['bin2_end'], row['bin1_end']
            return row

        return df.apply(swap_if_needed, axis=1)

    from pyliftover import LiftOver

    # Initialize LiftOver object for mm9 to mm10
    lo = LiftOver('mm9', 'mm10')

    jesse_ye_treg_loops = jesse_ye_loopdata[(jesse_ye_loopdata['logFC']>0) & (jesse_ye_loopdata['fdr'] < .05)]
    jesse_ye_tcon_loops = jesse_ye_loopdata[(jesse_ye_loopdata['logFC']<0) & (jesse_ye_loopdata['fdr'] < .05)]

    jesse_ye_treg_loops = liftover_df(jesse_ye_treg_loops, lo)
    jesse_ye_tcon_loops = liftover_df(jesse_ye_tcon_loops, lo)

    jesse_ye_treg_loops = swap_bins(jesse_ye_treg_loops)
    jesse_ye_tcon_loops = swap_bins(jesse_ye_tcon_loops)

    our_treg_loops = loopdict['Treg'].to_dataframe(header=None).iloc[:, :6]
    our_treg_loops.columns = ['chromosome', 'bin1_start', 'bin1_end', 'chrom2', 'bin2_start', 'bin2_end']
    our_treg_loops['chromosome'] = 'chr' + our_treg_loops['chromosome']

    our_tcon_loops = loopdict['Tcon'].to_dataframe(header=None).iloc[:, :6]
    our_tcon_loops.columns = ['chromosome', 'bin1_start', 'bin1_end', 'chrom2', 'bin2_start', 'bin2_end']
    our_tcon_loops['chromosome'] = 'chr' + our_tcon_loops['chromosome']

    jesse_ye_treg_loops['type'] = 'jesse_ye'
    our_treg_loops['type'] = 'ours'


    jesse_ye_tcon_loops['type'] = 'jesse_ye'
    our_tcon_loops['type'] = 'ours'

    # final_loop_df = pd.concat([jesse_ye_treg_loops, our_treg_loops], axis=0).reset_index(drop=True)

    comparisons = {
        'Tcon-specific loops' : [jesse_ye_tcon_loops, our_tcon_loops],
        'Treg-specific loops' : [jesse_ye_treg_loops, our_treg_loops],
        # 'JY Treg-Our Tcon' : [jesse_ye_treg_loops, our_tcon_loops],
        # 'JY Tcon-Our Treg' : [jesse_ye_tcon_loops, our_treg_loops],
    }
    final_overlap_results = []
    for comparison, (jy_loops, our_loops) in comparisons.items():
        final_loop_df = pd.concat([jy_loops, our_loops], axis=0).reset_index(drop=True)
        
        # Calculate delta_bin1 and delta_bin2
        delta_bin1 = np.abs(np.add.outer(final_loop_df['bin1_start'].values, -final_loop_df['bin1_start'].values))
        delta_bin2 = np.abs(np.add.outer(final_loop_df['bin2_start'].values, -final_loop_df['bin2_start'].values))
        same_chrom = np.equal.outer(final_loop_df['chromosome'].values, (final_loop_df['chromosome']).values)
        
        # Calculate overlap
        overlap_matrix = (same_chrom & (delta_bin1 < 25_005*3) & (delta_bin2 < 25_005*3))
        
        # Create the graph
        G = nx.Graph()
        
        # Add nodes with attributes
        for idx, row in final_loop_df.iterrows():
            G.add_node(idx, type=row['type'])
        
        # Add edges according to the overlap matrix
        for i in range(len(final_loop_df)):
            for j in range(i+1, len(final_loop_df)):
                if overlap_matrix[i, j]:
                    G.add_edge(i, j)
        
        # Find connected components
        connected_components = list(nx.connected_components(G))
        
        # Calculate fractions
        vals = []
        for c, component in enumerate(connected_components):
            for node in component:
                vals.append([c, G.nodes[node]['type']])
        
        vals = pd.DataFrame(vals, columns=['component', 'type'])
        unstack_df = (vals.value_counts(['component', 'type']).unstack().fillna(0))
        
        frac_jesse_with_ours = 1 - (unstack_df['jesse_ye'][(unstack_df['jesse_ye']>0) & (unstack_df['ours']>0)].sum())/(unstack_df['jesse_ye'].sum())
        frac_ours_with_jesse = 1 - (unstack_df['ours'][(unstack_df['jesse_ye']>0) & (unstack_df['ours']>0)].sum())/(unstack_df['ours'].sum())

        final_overlap_results.append([comparison, frac_jesse_with_ours, frac_ours_with_jesse])

    data = pd.DataFrame(final_overlap_results, columns = ['comparison', 'Loops in Liu et al. \n2023 not our loopset',
                                                'Loops in our loopset not in\n Liu et al. 2023 '])
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    sns.barplot(data=data.melt('comparison'), x='comparison', y='value', hue='variable',
            palette=['lightgray', 'tab:blue'])
    plt.xlabel("Loop set")
    plt.ylabel("Fraction")
    plt.title("Fraction unique loops (not in other loopset)")
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
    plt.xticks(rotation=20)
    fig.savefig('./plots/paper/s2/frac_unique_loops.pdf', bbox_inches='tight')    


def process_deseq2_output():
    loops = []
    pco = .05
    for file in glob.glob('./final_loops/DESEQ/output/all_merged_loops_all_bins_Treg_vs_Tcon_thresh*'):
        thresh = file.split("thresh=")[-1].split(".csv")[0]
        df = pd.read_csv(file, sep=' ')
        ns_loops = df[df['padj'] > pco].copy()
        loops_up = df[(df['padj'] < pco) & (df['log2FoldChange'] > 0)].copy()
        loops_down = df[(df['padj'] < pco) & (df['log2FoldChange'] < 0)].copy()
        all_loops = pd.concat([ns_loops, loops_up, loops_down], axis=0)
        
        ns_loops, loops_up, loops_down, all_loops = map(adjust_index_to_be_5kb, [ns_loops, loops_up, loops_down, all_loops])
        cond1, cond2 = file.split("_vs_")
        cond1 = cond1.split("bins_")[1]
        cond2 = cond2.split("_thresh")[0]
        
        all_loops = pd.concat([get_ancs_from_loops(all_loops), all_loops], axis=1)
        ns_loops = pd.concat([get_ancs_from_loops(ns_loops), ns_loops], axis=1)
        loops_up = pd.concat([get_ancs_from_loops(loops_up), loops_up], axis=1)
        loops_down = pd.concat([get_ancs_from_loops(loops_down), loops_down], axis=1)
        
        os.makedirs(f'final_loops/processed_DESEQ/thresh={thresh}', exist_ok=True)
        all_loops.to_csv(f'./final_loops/processed_DESEQ/thresh={thresh}/all_loops.csv', sep='\t', index=None, header=None)
        ns_loops.to_csv(f'./final_loops/processed_DESEQ/thresh={thresh}/ns_loops.csv', sep='\t', index=None, header=None)
        loops_up.to_csv(f'./final_loops/processed_DESEQ/thresh={thresh}/{cond1}_loops.csv', sep='\t', index=None, header=None)
        loops_down.to_csv(f'./final_loops/processed_DESEQ/thresh={thresh}/{cond2}_loops.csv', sep='\t', index=None, header=None)
    
        subdf_dict = {
            "ns": ns_loops,
            cond1: loops_up,
            cond2: loops_down,
            'all': all_loops,
        }
        for cond, df in subdf_dict.items():
            anc_granges = loop_grange_list_to_anc_granges(df.index)
            pbt.BedTool(anc_granges).saveas(f'./final_loops/processed_DESEQ/thresh={thresh}/{cond}_ancs.csv')
    return None


import matplotlib.pyplot as plt
import pybedtools as pbt
from matplotlib_venn import venn3_unweighted
from aux_functions import *

def make_venn_diagram_plot_s3():
    threshes = [0, 1]
    figs = []
    for thresh in threshes:
        ns_ancs = pbt.BedTool(f'final_loops/processed_DESEQ/thresh={thresh}/ns_ancs.csv')
        Treg_ancs = pbt.BedTool(f'final_loops/processed_DESEQ/thresh={thresh}/Treg_ancs.csv')
        Tcon_ancs = pbt.BedTool(f'final_loops/processed_DESEQ/thresh={thresh}/Tcon_ancs.csv')
    
        ns_ancs = bedtool_to_grange_set(ns_ancs)
        Treg_ancs = bedtool_to_grange_set(Treg_ancs)
        Tcon_ancs = bedtool_to_grange_set(Tcon_ancs)
        
        fig, ax = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi=200)
        venn3_unweighted([ns_ancs, Treg_ancs, Tcon_ancs], ['NS', "Treg", "Tcon"])
        title = f'Loop anchors using LFC threshold = {thresh}'
        ax.set_title(title)
        fig.savefig(f'./plots/paper/s3/venn_diagram_thresh={thresh}.pdf', bbox_inches='tight', dpi=300)



import seaborn as sns
import pybedtools as pbt
import matplotlib.pyplot as plt

def make_loop_lfc_cdf_plot(resting_gene_bedtool):
    threshes = [0, 1]
    for thresh in threshes:
        ns_ancs = pbt.BedTool(f'final_loops/processed_DESEQ/thresh={thresh}/ns_ancs.csv')
        Treg_ancs = pbt.BedTool(f'final_loops/processed_DESEQ/thresh={thresh}/Treg_ancs.csv')
        Tcon_ancs = pbt.BedTool(f'final_loops/processed_DESEQ/thresh={thresh}/Tcon_ancs.csv')
        
        ns_lfcs =   -get_col(resting_gene_bedtool.intersect(ns_ancs, u=True), 6).astype(float)
        Treg_lfcs = -get_col(resting_gene_bedtool.intersect(Treg_ancs, u=True), 6).astype(float)
        Tcon_lfcs = -get_col(resting_gene_bedtool.intersect(Tcon_ancs, u=True), 6).astype(float)
        
        fig, ax = plt.subplots(1, 1, figsize=(40*mm, 40*mm), dpi=200)
        sns.ecdfplot(ns_lfcs, color='black', label=f'Genes in NS anchors, n={len(ns_lfcs)}', linewidth=.5, ax=ax)
        sns.ecdfplot(Treg_lfcs, color='red', label=f'Genes in Treg anchors, n={len(Treg_lfcs)}', linewidth=.5, ax=ax)
        sns.ecdfplot(Tcon_lfcs, color='blue', label=f'Genes in Tcon anchors, n={len(Tcon_lfcs)}', linewidth=.5, ax=ax)
        
        ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=False)
        title = f"RNA-seq LFC (Loop |LFC| > {thresh})"
        ax.set_title(title)
        ax.set_xlim([-6, 4])
        ax.set_xlabel("RNA-seq LFC (Treg ÷ Tcon)")
        ax.set_ylabel(f"CDF")
        p_treg = format_pvalue(scipy.stats.ranksums(ns_lfcs, Treg_lfcs)[1])
        p_tcon = format_pvalue(scipy.stats.ranksums(ns_lfcs, Tcon_lfcs)[1])
    
        plt.text(1, 0, f'p<{p_treg}', va='bottom', ha='right', transform=ax.transAxes, color='red', fontsize=6)
        plt.text(1, .08, f'p<{p_tcon}', va='bottom', ha='right', transform=ax.transAxes, color='blue', fontsize=6)

        fig.savefig(f'./plots/paper/s3/rna_diagram_thresh={thresh}.pdf', bbox_inches='tight', dpi=300)
        
def plot_apa(full_pileup_dict, idx=None, res=5000):
    res = res / 1e3
    fig, axs = init_subplots_exact(6, 2, fgsz=(15*mm, 15*mm), dpi=100, sharey=True,
                                   xspace=1.4, yspace=1.4)
    for c, key in enumerate(full_pileup_dict):
        ax = axs[c]
        plt.sca(axs[c])
        if key == 'OGYuri.Tn.3':
            key = 'OGYuri.Tn.2'
        if idx is None:
            mats = full_pileup_dict[key]
        else:
            mats = full_pileup_dict[key][idx]
        mat = np.log2(np.nanmean(mats, axis=0))
        n = len(mat)//2
        m = len(mats)
        colorbar = axs[c].matshow(mat,             
                       extent = [-res*n, res*n, -res*n, res*n],
                       cmap='gist_heat_r', vmin=0, vmax=6)
        ax.set_xticks([-res*n, 0, res*n])
        ax.set_yticks([-res*n, 0, res*n])
        ax.tick_params(bottom=True, labelbottom=True, top=False, labeltop=False, left=True, labelleft=True)
        if c % 3 != 0:
            ax.tick_params(left=False, labelleft=False)
        if c < 3:
            ax.tick_params(bottom=False, labelbottom=False)

        ax.text(-res*n, -res*(n-1), f'n={m}', fontsize=6)
        ax.set_title(key.replace("OGYuri.", '').replace('.', ' rep'))
        ax.grid(False)
        
        if c % 3 == 2:
            cbar_ax = ax.inset_axes(transform = axs[c].transAxes, bounds = (1.1, 0, .05, 1))
            plt.colorbar(colorbar, cax=cbar_ax, shrink=.5)
            cbar_ax.grid(False)
            fig.add_axes(cbar_ax)
            plt.sca(cbar_ax)
            plt.yticks(fontsize=4)  
    return fig

import pandas as pd
import numpy as np
import pybedtools as pbt
import cooler
import matplotlib.pyplot as plt
from apa_figures import apa_diag_nonorm
from aux_functions import add_chr_to_bedtool

def make_tad_plot2(all_boundaries):

    cooldict = {
        'Tn.1' : cooler.Cooler('/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_cools/OGYuri.Tn.1.cool::/resolutions/5000'),
        'Tn.2' : cooler.Cooler('/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_cools/OGYuri.Tn.2.cool::/resolutions/5000'),
        'Tn.3' : cooler.Cooler('/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_cools/OGYuri.Tn.3.cool::/resolutions/5000'),
        'Treg.1' : cooler.Cooler('/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_cools/OGYuri.Treg.1.cool::/resolutions/5000'),
        'Treg.2' : cooler.Cooler('/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_cools/OGYuri.Treg.2.cool::/resolutions/5000'),
        'Treg.3' : cooler.Cooler('/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/zoomified_cools/OGYuri.Treg.3.cool::/resolutions/5000'),

    }

    apa_results, places = apa_diag_nonorm(all_boundaries, cooldict, wsz=40)
    shift = int(40 * 5_000 // 1e3)
    fig, axs = init_subplots_exact(6, 2, fgsz=(20*mm, 20*mm), dpi = 100, sharey=True)
    for c, key in enumerate(apa_results.keys()):
        plt.sca(axs[c])
        ax = axs[c]
        n = len(apa_results[key])
        mat = np.nanmean(apa_results[key], axis=0)
        mat[np.isnan(mat)] = 0
        mat += mat.T
        axs[c].matshow(np.log2(mat), cmap='gist_heat_r',
                    extent=[-shift, shift, shift, -shift]
                    )
        plt.grid(False)
        plt.xticks([-shift, 0, shift])
        plt.yticks([-shift, 0, shift])
        ax.set_xticklabels([-shift, 'Boundary', shift])
        ax.set_yticklabels([-shift, 'Boundary', shift])
        ax.tick_params(bottom=True, top=False, labeltop=False, labelbottom=True)
        if c % 3 != 0:
            ax.tick_params(labelleft=False, left=True)
        if c < 3:
            ax.tick_params(labelbottom=False, bottom=True)
        plt.title(key.replace(".", " rep").replace("Tn", "Tcon"))
        plt.text(0.02, 0.02, s=f'n={n}', transform=ax.transAxes, fontsize=6)
        fig.savefig('./plots/paper/s3/tad_pileup_replicates.pdf', bbox_inches='tight')

def make_looplength_plot(loopdict, loop_colordict):
    fig, ax = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi = 300)
    for i, loops in loopdict.items():
        l1 = get_col(loops, 1).astype(int)
        l2 = get_col(loops, 4).astype(int)
        sns.histplot(l2 - l1, ax = ax, color=loop_colordict[i], label=i, fill=False, element='step', stat='density', linewidth=1,
                    )
    ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))
    plt.legend(bbox_to_anchor=(1, 1))
    plt.xlabel("Distance")
    plt.title("Loop length")
    fig.savefig('./plots/paper/s3/loop_length.pdf',  bbox_inches='tight', dpi=300)


import pandas as pd
import pybedtools as pbt
import matplotlib.pyplot as plt

def make_tss_lineplot(resting_gene_lfcs, anchordict):
    transcribed_genes = resting_gene_lfcs[resting_gene_lfcs['itemRgb'] > 50]
    v = transcribed_genes['itemRgb']
    quantile_cut = pd.qcut(v, np.linspace(0, 1, 10))
    
    # Accessing the codes (integer representation) of the bins
    categorical_data = quantile_cut.cat
    codes = categorical_data.codes
    
    ls = []
    readcount = []
    quantiles = []
    for i in sorted(codes.unique()):
        subdf = transcribed_genes.loc[codes == i].reset_index()
        tmp_gene_bedtool = pbt.BedTool.from_dataframe(subdf[['chrom', 'start', 'end']])
        ls.append(bedprop(tmp_gene_bedtool, anchordict['All']))
        readcount.append(subdf['itemRgb'].mean())
        quantiles.append(i)
    
    fig, ax = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi=300)
    plt.plot(arr(quantiles) / np.max(arr(quantiles)), ls, marker='o', markersize=2)
    plt.xlabel("TSS by gene expression")
    plt.ylabel("Fraction")
    plt.title("Frac. TSS w/ Loop")
    plt.yticks([0, .25, .5, .75, 1])
    plt.xticks([0, .25, .5, .75, 1])
    ax.set_xticklabels(['Least expr.', '', '', '', 'Most expr.'])
    fig.savefig('./plots/paper/s3/frac_genes_in_loops.pdf', bbox_inches='tight', dpi=300)
    plt.show()

def loop_tad_overlap(loopdict, tad_boundaries):
    data = pd.DataFrame()
    for key in ['NS', 'Treg', 'Tcon']:
        loops = loopdict[key]
        loop_tad_overlap = loops.pair_to_bed(remove_chr_bedtool(pbt.BedTool.from_dataframe(tad_boundaries)), type='ispan').to_dataframe(
            header=None
        )
        values = list(loop_tad_overlap.iloc[:, :6].value_counts().reset_index()['count'])
        values = values + [0] * (len(loops)-len(values))
        values = pd.Series(values)
        values[values > 1] = 2
        
        values = pd.Series(values)
        values[values > 1] = 2
        values = values.value_counts(normalize=True).loc[[0, 1, 2]]
        data[key + " loops"] = values


    for key in ['NS', 'Treg', 'Tcon']:
        loops = loopdict[key]
        loop_tad_overlap = loops.pair_to_bed(remove_chr_bedtool(pbt.BedTool.from_dataframe(tad_boundaries
                                                                                        ).shift(s=200_000, genome='mm10')), type='ispan').to_dataframe(
            header=None
        )
        values = list(loop_tad_overlap.iloc[:, :6].value_counts().reset_index()['count'])
        values = values + [0] * (len(loops)-len(values))
        values = pd.Series(values)
        values[values > 1] = 2
        
        values = pd.Series(values)
        values[values > 1] = 2
        values = values.value_counts(normalize=True).loc[[0, 1, 2]]
        data[key + " loops (200kb shift)"] = values

    data.index = [0, 1, '>1']
    fig, axs = init_subplots_exact(2, 1, fgsz=(30*mm, 30*mm), dpi = 150)

    data.loc[:, ~data.columns.str.contains("shift")].T.plot.bar(stacked=True, color=['tab:blue', 'lightblue', 'lightgray'], ax=axs[0])
    data.loc[:, data.columns.str.contains("shift")].T.plot.bar(stacked=True, color=['tab:blue', 'lightblue', 'lightgray'], ax=axs[1])
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left', title='# TADs  crossed')
    fig.savefig('./plots/paper/s3/loop_tad_overlap.pdf', bbox_inches='tight')

import numpy as np
import matplotlib.pyplot as plt
from plotting_functions import add_yaxis_labels, init_subplots_exact

def make_chip_compartment_plots(my_treg_comp_50kb, bw_val_df_all_50kb):
    fig, axs = init_subplots_exact(5, 1, fgsz=(40*mm, 40*mm), dpi=200, sharey=True, xspace=1.5)
    for c, key in enumerate(['Treg H3K27ac', 'Treg H3K27me3', 'CD4 H3K9me3', 'Treg H3K4me3', 'Treg H3K4me1']):
        plt.sca(axs[c])
        name = key.split(" ")[1]
        plt.scatter(
            np.ravel(my_treg_comp_50kb)[:49264],
            np.ravel(np.log2(bw_val_df_all_50kb[key]) + 1)[:49264],
            s=.05, linewidth=0,
            rasterized=True
        )
        add_yaxis_labels("Low", "High", plt.gca(), fontsize=8, x=-.15)
        add_xaxis_labels("B", "A", plt.gca(), fontsize=8, )
        plt.ylabel(f"{name}")
        plt.xlabel("Compartment")
        plt.title(f"Compartment vs. {name}")
    return fig


import seaborn as sns
import matplotlib.pyplot as plt


def make_compartment_jointplot(my_treg_comp_50kb, my_tcon_comp_50kb):
    r, _ = samesize_nonan_test(my_treg_comp_50kb, my_tcon_comp_50kb)
    f = sns.jointplot(x=my_treg_comp_50kb, y=my_tcon_comp_50kb, s=.1, height=40*mm,
                      rasterized=True)
    f.ax_marg_x.grid(False)
    f.ax_marg_y.grid(False)
    f.set_axis_labels("Treg", "Tcon")
    f.ax_joint.set_xticks([-1.25, 0, 1.25])
    f.ax_joint.set_yticks([-1.25, 0, 1.25])
    f.ax_joint.set_xticklabels(['B', 0, 'A'])
    f.ax_joint.set_yticklabels(['B', 0, 'A'])
    f.fig.dpi = 300
    f.ax_joint.set_axisbelow(True)
    f.fig.suptitle("Compartment Score", y=1.2)
        
    # Annotate with Pearson R
    f.ax_joint.annotate(f'Pearson R = {r:.2f}', xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='center',
    fontsize=6)

    # Plot one point in red
    # Assuming you want to plot the first point in red, adjust the index as needed
    # f.ax_joint.scatter(my_treg_comp_50kb[278*5+7-10], my_tcon_comp_50kb[278*5+7-10], s=.1, color='red')

    
    f.savefig('plots/paper/s6/compartment_jointplot.pdf', bbox_inches='tight', dpi=1200)
    plt.show()



import numpy as np
import matplotlib.pyplot as plt
from supp_figures_plotting_functions import init_subplots_exact

def balanced_supplement_plots(deseq_lfc_mat, deseq_effect_mat):
    fig, ax = init_subplots_exact(2, 1, fgsz=(40*mm, 40*mm), dpi=200, xspace=1.4)
    
    xs, ys = [], [] 
    for i in range(len(deseq_lfc_mat)):
        vs = np.diag(deseq_lfc_mat, k=i)
        y = np.nanmedian(vs)
        xs.append(i)
        ys.append(y)
    ax[0].scatter(xs, ys, s=4)
    ax[0].set_title("Median, LFC")
    ax[0].set_xlabel("Diag (250kb)")
    ax[0].set_ylabel("Median")
    
    xs, ys = [], [] 
    for i in range(len(deseq_effect_mat)):
        vs = np.diag(deseq_effect_mat, k=i)
        y = np.nanmedian(vs)
        xs.append(i)
        ys.append(y)
    
    xs, ys = np.array(xs), np.array(ys)
    bad = ys == 0
    ax[1].scatter(xs[~bad], ys[~bad], s=4)
    ax[1].set_title("Median, Wald")
    ax[1].set_xlabel("Diag (250kb)")
    ax[1].set_ylabel("Median")
    
    for a in ax:
        a.set_axisbelow(True)
    
    fig.savefig('./plots/paper/s7/balanced_wald.pdf', bbox_inches='tight')
    plt.show()


import numpy as np
import matplotlib.pyplot as plt
import math

def deseq_lfc_vs_compartment_fc(all_ind_to_region, chrom_to_start, chrom_to_end, my_treg_comp, deseq_lfc_mat, deseq_pval_mat, ind_to_gene):
    index_list = [161, 243, 278, 260, 3208, 828, 3207, 354, 4379, 4380, 2732, 2035, 6501, 5699, 5700]
    n = len(index_list)
    cols = 3
    rows = math.ceil(n / 3)
    fig, axs = init_subplots_exact(n, cols, fgsz=(30*mm, 30*mm), dpi=200, xspace=1.4, yspace=1.6)
    axs = np.ravel(axs)
    for c, i in enumerate(index_list):
        chrom = all_ind_to_region[i][0]
        s, e = chrom_to_start[chrom], chrom_to_end[chrom]
        ax = axs[c]
        x = my_treg_comp[s:e]
        y = deseq_lfc_mat[i, s:e]
        p = deseq_pval_mat[i, s:e]
        ax.scatter(x, y, zorder=3, s=4, c=(p < .05), cmap='coolwarm', vmin=-1, vmax=1)
        r = samesize_nonan_test(x, y)[0]
        ax.set_title(f'Gene={get_name(i, ind_to_gene)};\nr={r.round(2)};')
        ax.set_xlabel("Treg comp")
        ax.set_ylabel("DESeq2 LFC")
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-5.5, 5.5])
    plt.tight_layout()
    fig.savefig('./plots/paper/s7/compartment_lfc_corr.pdf', bbox_inches='tight')
    plt.show()


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def deseq_lfc_vs_ab_compartment_correlations(all_ind_to_region, chrom_to_start, chrom_to_end, my_treg_comp, deseq_lfc_mat, deseq_effect_mat):
    rs = []
    for i in range(len(all_ind_to_region)):
        chrom = all_ind_to_region[i][0]
        s, e = chrom_to_start[chrom], chrom_to_end[chrom]
        x = my_treg_comp[s:e]
        y = deseq_lfc_mat[i, s:e]
        r = samesize_nonan_test(x, y)[0]
        rs.append(r)
    rs = arr(rs)
    
    goodinds = np.nansum(np.abs(deseq_effect_mat) > 3, axis=1) > 10
    
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi=200, as_list=True)
    ax = axs[0]
    sns.histplot(rs[goodinds], ax=ax, zorder=2, alpha=.9)
    ax.set_xlabel("Pearson Rs")
    ax.set_title("Pearson R: \nDESeq LFCs vs. A/B Compartment")
    fig.savefig('./plots/paper/s7/lfc_compartment_r.pdf', bbox_inches='tight')
    plt.show()


def n_genes_per_bin(ind_to_gene, all_ind_to_region):
    lens = []
    l = np.arange(len(all_ind_to_region))
    for i in l:
        n = ind_to_gene.get(i, [])
        lens.append(len(n))
    
    lens = arr(lens)
    lens[lens > 20] = 20
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi=200)
    sns.histplot(lens)
    plt.xlim([-1, 21])
    plt.xlabel("# Genes per 250kb bin")
    plt.title("Genes per 250kb bin")
    fig.savefig('./plots/paper/s7/genes_per_250kb_bin.pdf', bbox_inches='tight')    



import numpy as np
import matplotlib.pyplot as plt
import statsmodels
import statsmodels.nonparametric.smoothers_lowess
from aux_functions import *

def ikzf2_histone_correlations(bw_val_df_all_250kb, deseq_effect_mat, deseq_pval_mat):
    n = len(bw_val_df_all_250kb)
    xlims = [[0, 3.5], [0, 3.5], [0, 1.5], [0, 3.5], [0, 1.5], [0, 3.5]]
    fig, axs = init_subplots_exact(6, 2, dpi=200, fgsz=(40*mm, 40*mm), space=1.4, sharey=True)
    
    for c, key in enumerate(['Treg H3K27ac', 'Tcon H3K27ac', 'CD4 H3K9me3', 'Treg H3K4me1', 'Treg H3K4me3', 'Treg H3K27me3']):
        vector = np.ravel(bw_val_df_all_250kb[key])
        name = key.split(" ")[1]
        plt.sca(axs[c])
        x = vector[:782]
        y = deseq_effect_mat[278, :782]
        x = np.ravel(x)
        
        p = deseq_pval_mat[278, :782]
        bad = np.isnan(x) | np.isnan(y)
        x, y = x[~bad], y[~bad]
        p = p[~bad]
        x[x > 3] = 3
    
        r = samesize_nonan_test(x, y)[0]
        
        smooth_x, smooth_y = statsmodels.nonparametric.smoothers_lowess.lowess(x, y).T
        plt.scatter(x, y, zorder=3, s=4, c=(p < .05), cmap='coolwarm', vmin=-1, vmax=1)
        plt.plot(smooth_x, smooth_y, color='black', zorder=3, linewidth=.5, linestyle='--')
        plt.title(name)
        plt.ylabel("Ikzf2 DESeq2 LFC")
        plt.xlabel(f"{name}")
        plt.title(key + f'\nr={r.round(2)}')
        add_xaxis_labels("Low", "High", plt.gca(), fontsize=6)
        add_yaxis_labels("Tcon", "Treg", plt.gca(), fontsize=6, x=-.15)
        plt.xlim(xlims[c])
    
    fig.savefig('./plotspaper/s7/ikzf2_histone_figure.pdf', bbox_inches='tight')
    plt.show()



def deseq2_reproducibility_heatmap(deseq_pval_mat, deseq_effect_mat):
    all_read_counts = pd.read_csv('/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/deseq_interactions/all_intra_chromosomal.exp=None.rebal=False.csv',
                              sep='\t')
    all_read_counts = all_read_counts[all_read_counts['ind'].str.contains('278')]
    all_read_counts = all_read_counts[all_read_counts['chrom']==1]

    
    p = deseq_pval_mat[all_read_counts['rows'], all_read_counts['cols']]
    effect = deseq_effect_mat[all_read_counts['rows'], all_read_counts['cols']]
    
    cols = sns.color_palette('bwr', n_colors=5)
    row_colors = []
    for x in (p < .05)*np.sign(effect):
        if x == -1:
            row_colors.append(cols[0])
        elif x == 1:
            row_colors.append(cols[-1])
        else:
            row_colors.append(cols[2])
    sign_colors = pd.Series(row_colors)
        
    # cols = sns.color_palette('bwr', as_cmap=True)
    # row_colors = []
    # for x in effect:
    #     x = np.clip(x, -7.99, 7.99)
    #     x = (x+8)/16
    #     if np.isnan(x):
    #         row_colors.append(np.nan)
    #     else:
    #         row_colors.append(cols(x))
    # effect_colors = pd.Series(row_colors)
    
    
    # cols = sns.color_palette('gray_r', as_cmap=True)
    # row_colors = []
    # for x in p:
    #     x = -np.log10(x)/2
    #     if np.isnan(x):
    #         x = 0
    #     else:
    #         x = min(.999, x)
    #     row_colors.append(cols(x))
    # p_colors = pd.Series(row_colors)
    
    
    row_colors = pd.concat([sign_colors], axis=1)
    row_colors.columns = ['Differential']
    data = all_read_counts.iloc[:, :6].apply(zscore, axis=1)
    idx = ~(data.isna().any(axis=1) | row_colors.isna().any(axis=1))
    row_colors.index = data.index
    
    g = sns.clustermap(data[idx], cmap='coolwarm', vmin=-2, vmax=2, row_colors=row_colors[idx],
                      yticklabels=False)
    return g


### Figure S7

def make_chromosomewide_deseq_plot(deseq_effect_mat, deseq_pval_mat, chrom, 
                                   chrom_to_start, chrom_to_end, inds_to_add = [], dpi = 300):
    s, e = chrom_to_start[chrom], chrom_to_end[chrom]
    fig, ax = init_subplots_exact(1, 1, fgsz=(60*mm, 60*mm), dpi=dpi)
    pos = ax.matshow(deseq_effect_mat[s:e, s:e]*(deseq_pval_mat[s:e, s:e] < .05), cmap=cm.bwr, vmin=-8, vmax=8,
                    interpolation='none')
    if chrom == '1':
        ax.set_xticks([0, 200, 277, 400, 600, 781])
        ax.set_yticks([0, 200, 277, 400, 600, 781])
    
        ticks = list(ax.get_xticks()*250_000//1_000_000)
        ticks[2] = "Ikzf2"; 
    else:
        ticks = list(np.arange(0, e-s, 200))
        for ind, name in inds_to_add:
            d = (ind-1-s)
            ticks.append(d)
        ticks = sorted(ticks)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ticks = list(ax.get_xticks()*250_000//1_000_000)
        for ind, name in inds_to_add:
            d = str(int((ind-1-s)*250_000//1e6))
            idx = np.where(arr(ticks).astype(str)==d)[0][0]
            ticks[idx] = name
            
    ax.tick_params(axis='x', bottom=True, top=False, labeltop=False, labelbottom=True)
    ax.tick_params(axis='y', left=False, right=True, labelright=True, labelleft=False)

    ax.set_xticklabels(ticks)
    ax.set_yticklabels(ticks)
    
    ax.set_xlabel(f"chr{chrom} (Mb)",)
    ax.set_ylabel(f"chr{chrom} (Mb)",)
    ax.yaxis.set_label_position("right")
    ax.set_title(f"Treg Hi-C ÷ Tcon Hi-C, chr{chrom}")
    
    fig.canvas.draw()
    
    newax = ax.inset_axes(transform=ax.transAxes,
            bounds = (-.2, .25, .03, .5))
    plt.colorbar(pos, cax=newax)
    newax.grid(False)
    newax.text(0.5, -.1, s='Tcon', ha = 'center', va = 'center', transform=newax.transAxes)
    newax.text(0.5, 1.1, s='Treg', ha = 'center', va = 'center', transform=newax.transAxes)
    newax.tick_params(axis='y', left=True, right=False, labelright=False, labelleft=True)
    newax.set_yticks([-8, -4, 0, 4, 8])
    newax.set_ylabel("DESeq2 Wald Statistic")
    newax.yaxis.set_label_position("left")
    
    newax = ax.inset_axes(bounds=(0, -.22, 1, .1), 
                       transform=ax.transAxes,)
    newax.set_xticks([])
    inds = np.arange(0, e-s)
    newax.set_xlim(0, e-s)
    newax.set_ylim([-250, 250])

    treg_specific = np.sum((deseq_pval_mat[s:e, s:e] < .05) & (deseq_effect_mat[s:e, s:e] > 0), axis=1)
    tcon_specific = np.sum((deseq_pval_mat[s:e, s:e] < .05) & (deseq_effect_mat[s:e, s:e] < 0), axis=1)
    newax.fill_between(inds, [0]*len(treg_specific), treg_specific, label='# Treg-specific\ninteractions', alpha=1, color='red')
    newax.fill_between(inds, [0]*len(tcon_specific), -tcon_specific, label='# Tcon-specific\ninteractions', alpha=1, color='blue')
    # newax.axis('off') 
    newax.set_yticks([-250, 250])
    # plt.yticks(fontsize=12)
    for spine in newax.spines.values():
        spine.set_visible(False)
    newax.spines['left'].set_visible(True)
    newax.tick_params(top='off', bottom='off', left=True, labelleft=True,
                      labelbottom='off')
    newax.grid(False)
    # newax.yaxis.set_ticks_position('none') 
    newax.legend(bbox_to_anchor=(1.05, .52), loc='center left', )
    plt.grid(False)
    print(chrom)
    path = f'./plots/paper/s8/{chrom}_treg_tcon_hicmat.pdf'
    fig.savefig(path, bbox_inches='tight')
    return fig


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def make_manhattan_plot(deseq_pval_mat, all_ind_to_region, chrom_to_start, chrom_to_end, ind_to_gene, bin_cutoff=0, label_cutoff=100, dpi=200):
    v = (np.triu(np.abs(deseq_pval_mat < .05), k=bin_cutoff) 
         + np.tril(np.abs(deseq_pval_mat < .05), k=-bin_cutoff)).sum(axis=1)
    xticks = []
    fig, axs = init_subplots_exact(1, 1, fgsz=(100*mm, 30*mm), dpi=dpi)
    chroms = pd.unique(get_col(all_ind_to_region, 0))
    
    for c, chrom in enumerate(chroms):
        s = chrom_to_start[chrom]
        e = chrom_to_end[chrom]
        xs = np.arange(s, e)
        color = 'lightblue' if c % 2 == 0 else 'tab:blue'
        axs.scatter(xs, v[s:e], color=color, s=8, linewidth=0)
        xticks.append((s + e) / 2)
    
    axs.set_xticks(xticks)
    axs.set_xticklabels(add_chr_to_list(chroms))
    plt.xticks(rotation=90)
    
    seen = set()
    o = np.argsort(-v)
    for i in o:
        if v[i] < label_cutoff:
            continue
        name = get_name(i, ind_to_gene)
        if name == 'None':
            continue
        if (np.abs((np.asarray(list(seen)) - i)) < 15).any():
            continue
        seen.add(i)
        axs.text(i, v[i], name, fontsize=6)
        axs.scatter(i, v[i], edgecolor='black', s=8, facecolor='none')
    
    axs.set_axisbelow(True)
    axs.set_xlabel("")
    axs.set_ylabel("# Differential Interactions")
    
    return fig



### Fig S10.5
from plotting_functions import plot_from_multiple_coolers_and_one_grange
def make_replicate_plot(replicate_cooldict, all_ind_to_region, ind_to_gene):
    inds = [[5145, 5216], [243, 278], [426, 772], [2179, 2341]]
    for ind1, ind2 in inds:
        name1, name2 = get_name(ind1, ind_to_gene), get_name(ind2, ind_to_gene)
        print(name1, name2)
        if ind2 == 772:
            name2 = 'Irf6'
        if ind2 == 2341:
            name2 = 'Tgfbr1'
        if ind1 == 5145:
            name1 = 'Gpr83/Izumo1r'
        fig = plot_from_multiple_coolers_and_one_grange(replicate_cooldict,
                                                tuple_to_grange(*all_ind_to_region[ind1]) + "|" + tuple_to_grange(*all_ind_to_region[ind2]),
                                                250_000*2, useSigma=False,
                                                        ylabel=name1,
                                                        xlabel=name2,
                                                );
        name1 = name1.replace("/", "_")
        fig.savefig(f'./plots/paper/s10/replicate_plot_{name1}_{name2}.pdf', bbox_inches='tight')



def generate_metadomain_venn_diagrams(deseq_pval_mat, all_intra_treg_metadomains, all_intra_tcon_metadomains, pcos = [.05, 1]):
    for c, pco in enumerate(pcos):
        if c == 0:
            idx = (deseq_pval_mat < pco)
        else:
            idx = (deseq_pval_mat < pco) & (deseq_pval_mat > pcos[c-1])
        
        treg_mega_df = get_metadomains_by_distance(all_intra_treg_metadomains * idx)
        tcon_mega_df = get_metadomains_by_distance(all_intra_tcon_metadomains * idx)
        
        print(treg_mega_df.shape)

        
        treg_metadomain_set = ind_mega_df_to_set(treg_mega_df)
        tcon_metadomain_set = ind_mega_df_to_set(tcon_mega_df)
        
        fig, axs = init_subplots_exact(1, 1, fgsz=(30*mm, 30*mm), dpi=200, sharex=True, space=1.5)
        my_venn2([tcon_metadomain_set, treg_metadomain_set], ['Tcon', 'Treg'], ax=axs)
        jacc_coeff = np.round(len(tcon_metadomain_set.intersection(treg_metadomain_set)) / len(tcon_metadomain_set.union(treg_metadomain_set)),
                              2)
        axs.set_title(f"Jaccard: {jacc_coeff}")
        fig.savefig(f'./plots/paper/s11/metadomain_venn_diagram_pco={pco}.pdf', bbox_inches='tight', dpi=300)
        plt.show()


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from volcano_plot import volcano_plot
from aux_functions import get_name

def generate_metadomain_volcano_plot(all_intra_metadomains, deseq_lfc_mat, deseq_pval_mat, ind_to_gene, ignore_set, output_path):
    n = len(all_intra_metadomains)
    X, Y = np.indices((n, n))
    indsoi = np.triu(all_intra_metadomains, k=8)
    inds_x, inds_y = X[indsoi], Y[indsoi]

    data = pd.DataFrame()
    data['lfc'] = deseq_lfc_mat[indsoi]
    data['pval'] = deseq_pval_mat[indsoi]
    data['genes'] = ['-'.join((get_name(x, ind_to_gene), get_name(y, ind_to_gene))) for x, y in zip(inds_x, inds_y)]

    fig, ax = plt.subplots(figsize=(70*mm, 70*mm), dpi = 150)
    volcano_plot(data['lfc'], data['pval'], data['genes'], ax, max_y=35, ylim=[0, 36], 
                 label_pval_cutoff=15, xlim=[-3, 3], ignore_set=ignore_set, fontsize=6,
                 pco = .05,
                 lfc_co=0,
                 pes=None)

    ax.set_xlabel(r'Log$_2$FC')
    ax.set_ylabel(r'-Log$_{10}$ FDR')
    ax.set_title(r'Megaloop Volcano Plot')
    fig.savefig(output_path, bbox_inches='tight')
    plt.show()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.gridspec as gridspec

def generate_metadomain_count_plot(deseq_effect_mat, all_intra_metadomains, cutoff=4):
    n_de_metadomains = np.sum((np.abs(deseq_effect_mat) > 4) & (all_intra_metadomains > 0), axis=1)
    n_treg_metadomains = np.sum((deseq_effect_mat > 4) & (all_intra_metadomains > 0), axis=1)
    n_tcon_metadomains = np.sum((deseq_effect_mat < -4) & (all_intra_metadomains > 0), axis=1)
    n_metadomains = np.sum((all_intra_metadomains > 0), axis=1)

    metadomain_count_df = pd.DataFrame([n_de_metadomains, n_metadomains, n_treg_metadomains, n_tcon_metadomains]).T
    metadomain_count_df.columns = ['# DE metadomains', '# metadomains', '# Treg metadomains', '# Tcon metadomains']

    col = '# metadomains'
    cutoff = 1000
    data = metadomain_count_df
    data = data[col].value_counts().reset_index()
    tmp = pd.DataFrame(np.asarray(metadomain_count_df['# metadomains']), columns=[col]).value_counts().reset_index()
    data = tmp
    data = data.set_index(col).reindex(np.arange(max(data[col]))).reset_index()
    data[col] = data[col].astype(int)
    # data.loc[cutoff, 'count'] += data.loc[cutoff+1:, 'count'].sum()
    # data = data.loc[:cutoff]
    # data = data.loc[1:]

    # Set up the gridspec to define the layout of the axes
    fig, axs = init_subplots_exact(2, 2, fgsz=(40*mm, 15*mm), y_ratios=[.5, 1], sharex=True)
    ax0, ax1 = axs

    # Plot the same data on both axes
    ax0.bar(data[col], data['count'], zorder=3)
    ax1.bar(data[col], data['count'], zorder=3)

    # Set the y-axis limits
    ax0.set_ylim(100, 5200)
    ax1.set_ylim(0, 100)
    ax1.set_yticks([0, 50, 100])
    # Hide the spines between ax and ax2
    ax0.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # ax0.xaxis.tick_top()
    ax0.tick_params(labeltop=False, bottom=False)  # don't put tick labels at the top

    # Set labels and title
    ax0.set_ylabel("Count")
    ax1.set_ylabel("Count")
    ax1.set_xlabel("# Megaloops")
    ax0.set_title("# Megaloops by bin")
    for ax in [ax0, ax1]:
        plt.sca(ax)
        plt.yticks()

    # Draw diagonal lines on the broken axis
    d = .015  # size of break mark
    kwargs = dict(transform=ax0.transAxes, color='k', clip_on=False)
    ax0.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax0.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    ax0.tick_params(labelbottom=False)
    ax0.set_ylabel("")
    fig.savefig('./plots/paper/s11/metadomains_by_bin.pdf', bbox_inches = 'tight')



def generate_metadomain_compartment_plot(all_intra_metadomains, my_treg_comp, ACOMPARTMENT_CUTOFF_LOOSE, ):
    tmp1, tmp2 = np.where(all_intra_metadomains)
    # fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 300)
    x = my_treg_comp[tmp1]
    y = my_treg_comp[tmp2]
    f = sns.jointplot(x = x, y = y, s=.1, zorder=3, height=50*mm,
                    rasterized=True);

    ax = f.ax_joint
    shade_x = np.linspace(ACOMPARTMENT_CUTOFF_LOOSE, 2, 100)
    shade_y = np.linspace(ACOMPARTMENT_CUTOFF_LOOSE, 2, 100)
    shade_X, shade_Y = np.meshgrid(shade_x, shade_y)

    # Create a mask for shading
    mask = (shade_X >= ACOMPARTMENT_CUTOFF_LOOSE) & (shade_Y >= ACOMPARTMENT_CUTOFF_LOOSE)

    # Apply the shading
    ax.contourf(shade_X, shade_Y, mask, levels=[0.5, 1], colors='lightgray', alpha=0.5, zorder=1, 
                label='A compartment (strict)')

    frac = ((x > ACOMPARTMENT_CUTOFF_LOOSE) & (y > ACOMPARTMENT_CUTOFF_LOOSE)).mean()
    ax.annotate(f'{frac:.2f}', xy=(0.95, 0.95), xycoords='axes fraction', ha='right', va='center',
                fontsize=6)
    
    frac = ((x > ACOMPARTMENT_CUTOFF_LOOSE) & (y < ACOMPARTMENT_CUTOFF_LOOSE)).mean()
    ax.annotate(f'{frac:.2f}', xy=(0.95, 0.05), xycoords='axes fraction', ha='right', va='center',
                fontsize=6)
    
    frac = ((x < ACOMPARTMENT_CUTOFF_LOOSE) & (y > ACOMPARTMENT_CUTOFF_LOOSE)).mean()
    ax.annotate(f'{frac:.2f}', xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='center',
                fontsize=6)

    frac = ((x < ACOMPARTMENT_CUTOFF_LOOSE) & (y < ACOMPARTMENT_CUTOFF_LOOSE)).mean()
    ax.annotate(f'{frac:.2f}', xy=(0.05, 0.05), xycoords='axes fraction', ha='left', va='center',
                fontsize=6)

    plt.axvline(ACOMPARTMENT_CUTOFF_LOOSE, color='black', linestyle='--')
    plt.axhline(ACOMPARTMENT_CUTOFF_LOOSE, color='black', linestyle='--')

    plt.xlabel("Anchor 1")
    plt.ylabel("Anchor 2")
    plt.xticks([-1, 0, 1])
    plt.yticks([-1, 0, 1])
    plt.gca().set_xticklabels(['B', '', 'A'])
    plt.gca().set_yticklabels(['B', '', 'A'])
    plt.title("Megaloop compartments", y=1)
    f.ax_marg_x.remove()
    f.ax_marg_x.grid(False)
    f.ax_marg_y.grid(False)

    f.fig.dpi=300
    return f.fig



import matplotlib
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from megaloop_functions import *

import matplotlib
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_chromosome_heatmap(chrom, all_intra_metadomains, gene_to_ind, SE_count, bw_val_df_250kb, my_treg_comp, chrom_to_start, chrom_to_end, cutoff=5,
                                dpi=150, n_clusters = 10):
    color1 = (190/255, 104/255, 197/255)
    color2 = (251/255, 251/255, 251/255)
    custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_cmap", [color2, color1])
    mm = 1 / 2.54 / 10

    genes_to_plot = {
        '1': {'Dusp10', 'Cxcr4', 'Cd247', 'Fasl', 'Platr22', 'Itpkb', 'Ly96', 'Il17f', 'Ikzf2',  'Ctla4', 'Idh1', 
              'Bcl2', 'Irf6', 'Lct', 'Cxcr4', 'Cdc73', 'Ptprc',  'Pou2f1', 'Dusp23', 'Tgfb2'}
    }

    all_intra_clustering = []


    intra_clustering = pd.DataFrame()
    s, e = chrom_to_start[chrom], chrom_to_end[chrom]
    size = e - s
    submat = pd.DataFrame(all_intra_metadomains[s:e, s:e] > 0)
    indsoi = submat.sum(axis=0) >= cutoff
    indsoi = np.where(indsoi)[0]
    submat = submat.loc[indsoi, indsoi]

    row_colors = [sns.color_palette('bwr', as_cmap=True)(x / size) for x in indsoi]
    o, clusters = make_order_and_cluster_optimal(submat, method='ward', metric='euclidean', n_clusters=n_clusters)
    clusters = rename_clusts_by_order(clusters, o)

    fig, axs = init_subplots_exact(2, 1, fgsz=(70*mm, 70*mm), dpi=dpi, space=1.4)
    sns.heatmap(submat.iloc[o, o], cmap=custom_cmap, vmin=0, vmax=1, cbar=False, ax=axs[1], rasterized=True)
    sns.heatmap(submat, cmap=custom_cmap, vmin=0, vmax=1, cbar=False, ax=axs[0], rasterized=True)

    genes_in_order = []
    for gene in genes_to_plot.get(chrom, []):
        inds = np.array(gene_to_ind.get(gene)) - s
        try:
            ind = inds[np.argmax(submat.sum(axis=1).loc[[x for x in inds if x in submat.index]])]
        except:
            continue
        matrix_ind = np.where(submat.index[o] == ind)[0][0]
        genes_in_order.append([gene, ind, matrix_ind])

    n = len(genes_in_order)
    if n > 0:
        genes_in_order = pd.DataFrame(genes_in_order).sort_values(2, ascending=False)
        for c, (_, row) in enumerate(genes_in_order.iterrows()):
            name = row[0]
            matrix_ind = row[2]
            textx, texty = 1.05, c / n
            axs[1].text(textx, texty + 1 / n / 2, name, transform=axs[1].transAxes, fontsize=8, va='center')
            xpos = 1
            ypos = matrix_ind / len(submat.index)
            axs[1].annotate('', xy=(textx, texty + 1 / n / 2), xycoords='axes fraction', xytext=(xpos, 1 - ypos),
                            arrowprops=dict(arrowstyle="-", lw=.5, color='red'))

    for a in axs:
        a.set_xticks([])
        a.set_yticks([])

    
    pos1 = axs[1].get_position()
    ax3 = fig.add_axes([pos1.x0 - .05, pos1.y0, 0.03, pos1.height])
    sns.heatmap((indsoi[o] / size)[:, None], cmap='YlGn', vmin=0, vmax=1, cbar=False, ax=ax3, rasterized=True)
    ax3.set_title('Position', rotation=40, fontsize=6, ha='left')
    ax3.set_yticks([])
    
    pos1 = axs[0].get_position()
    ax2 = fig.add_axes([pos1.x0 - .05, pos1.y0, 0.03, pos1.height])
    sns.heatmap((indsoi / size)[:, None], cmap='YlGn', vmin=0, vmax=1, cbar=False, ax=ax2, rasterized=True)
    ax2.set_yticks([])
    ax2.set_ylabel('Chromosome position')

    pos1 = axs[1].get_position()
    ax3 = fig.add_axes([pos1.x0 - .1, pos1.y0, 0.03, pos1.height])
    sns.heatmap(np.ravel(clusters[o])[:, None], cmap='tab10', cbar=False, ax=ax3, rasterized=False, vmin=0, vmax=10)
    ax3.set_title("Cluster", rotation=40, fontsize=6, ha='left')
    ax3.set_yticks([])

    pos1 = axs[1].get_position()
    ax3 = fig.add_axes([pos1.x0 - .15, pos1.y0, 0.03, pos1.height])
    sns.heatmap((SE_count[indsoi + s][o] > 0)[:, None], cmap='gray_r', vmin=0, vmax=1, cbar=False, ax=ax3, rasterized=True)
    ax3.set_title("SE", rotation=40, fontsize=6, ha='left')
    ax3.set_yticks([])

    pos1 = axs[1].get_position()
    ax3 = fig.add_axes([pos1.x0 - .2, pos1.y0, 0.03, pos1.height])
    sns.heatmap((np.log2(1 + bw_val_df_250kb['Treg H3K27me3']).iloc[indsoi + s].iloc[o]).values[:, None], cmap='bwr', vmin=-1, vmax=2, cbar=False, ax=ax3, rasterized=True)
    ax3.set_title("H3K27me3", rotation=40, fontsize=6, ha='left')
    ax3.set_yticks([])

    pos1 = axs[1].get_position()
    ax3 = fig.add_axes([pos1.x0 - .25, pos1.y0, 0.03, pos1.height])
    sns.heatmap((np.log2(1 + bw_val_df_250kb['Treg H3K27ac']).iloc[indsoi + s].iloc[o]).values[:, None], cmap='bwr', vmin=-1, vmax=2, cbar=False, ax=ax3, rasterized=True)
    ax3.set_title("H3K27ac", rotation=40, fontsize=6, ha='left')
    ax3.set_yticks([])

    pos1 = axs[1].get_position()
    ax3 = fig.add_axes([pos1.x0 - .3, pos1.y0, 0.03, pos1.height])
    sns.heatmap(my_treg_comp[indsoi + s][o][:, None], cmap='bwr', vmin=-1, vmax=1, cbar=False, ax=ax3, rasterized=True)
    ax3.set_title("Treg comp", rotation=40, fontsize=6, ha='left')
    ax3.set_yticks([])


    axs[0].set_title(f'Bins w/ ≥ {cutoff} metadomains (Chr{chrom})')
    axs[1].set_title(f'Bins w/ ≥ {cutoff} metadomains, reordered (Chr{chrom})')

    intra_clustering['ind'] = indsoi + s
    intra_clustering['chr'] = chrom
    intra_clustering['clust'] = clusters
    all_intra_clustering.append(intra_clustering)

    fig.savefig(f'./plots/paper/s12/heatmap_chr={chrom}_cutoff={cutoff}.pdf', bbox_inches='tight')
    fig.savefig(f'./plots/paper/s12/heatmap_chr={chrom}_cutoff={cutoff}.svg')
    return all_intra_clustering



import math
from matplotlib.ticker import FuncFormatter

def million_base_formatter(x, pos):
    return f'{int(x * 1e-6)} Mb'


# print(1)
def make_wrapped_plot(inter_and_intra_connections, cool, all_ind_to_region, chrom_to_start, 
                      chrom_to_end, og_ind=66, compare_with = '3', dpi = 100, ylabel='',
                      gene='Ikzf2', intra=True,):
    d = 20
    og_chrom = all_ind_to_region[og_ind][0]
    ind = (og_ind-chrom_to_start[og_chrom])*5
    m = cool.matrix().fetch(og_chrom, compare_with)
    res = cool.info['bin-size']
    salt, ealt = chrom_to_start[compare_with], chrom_to_end[compare_with]
    
    cols = m.shape[1]
    plotsize = 500
    n = math.ceil(cols/500)
    ratio = plotsize // (d*2)
    fig, axs = init_subplots_exact(n, n, fgsz=(ratio, 1), dpi=dpi)
    for i in range(n):
        cols,  cole = i*plotsize, i*plotsize+plotsize
        if cole > m.shape[1]:
            axs[i].remove()
            continue
        if compare_with == og_chrom:
            vmax = 3e-3
        else:
            vmax = 5e-4 

        ymin, ymax = (ind-20)*res, (ind+20)*res

        submat = m[ind-d:ind+d, cols:cole]
        axs[i].matshow(submat, cmap='gist_heat_r', vmax = vmax, aspect='equal',
                      extent = [cols*res, cole*res, ymax, ymin])
        axs[i].grid(False)
        for col in np.where(inter_and_intra_connections[og_ind, salt:ealt])[0]:
            x = (col+.5)*5
            if (x > cols) and (x < cole):
                axs[i].arrow(x*res, ymax-5*res, 0, -2*res, color='black', head_width=8*res, 
                             head_length=4*res, width=4*res)
        axs[i].set_xlim([cols*res, cole*res])
        formatter = FuncFormatter(million_base_formatter)
        axs[i].xaxis.set_major_formatter(formatter)
        axs[i].yaxis.set_major_formatter(formatter)

        axs[i].set_ylabel(ylabel)
        axs[i].set_ylabel(gene)
        # axs[i].set_yticklabels([])
        axs[i].set_yticks([ymax, (ymax + ymin)//2, ymin])
        if (cols < ind) & (cole > ind):
            x = (ind+2.5)*res
            plt.sca(axs[i])
            current_ticks = plt.xticks()[0]
            current_labels = [item.get_text() for item in plt.xticks()[1]]
            new_tick = x
            new_label = gene
            # Append the new tick and label
            if intra:
                updated_ticks = list(current_ticks) + [new_tick]
                updated_labels = current_labels + [new_label]
                plt.xticks(updated_ticks, updated_labels)
    return fig



from adjustText import adjust_text
def treg_tcon_inter_metadomains(all_inter_treg_metadomains, all_inter_tcon_metadomains, ind_to_gene,):
    treg_inter = all_inter_treg_metadomains.sum(axis=1)
    tcon_inter = all_inter_tcon_metadomains.sum(axis=1)
    x, y = treg_inter, tcon_inter

    fig, ax = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    scatter_with_pearson(x, y, s=8, ax=plt.gca(), zorder=3, color = 'skyblue', rasterized=True)

    text_labels = []

    genes_to_plot = set(list(np.where(y > 150)[0]) + list(np.where((x - y < -65))[0][::3]) + list(np.where(((x-y) > 20))[0]) + [9010])
    for i in genes_to_plot:
        gene = get_name(i, ind_to_gene)
        if gene not in ['Tox', 'Ikzf2', 'Ikzf1', 'Pik3ip1', 'Lef1', 'Tgfbr3', 'Runx1', 'Lef1', 'Lta']:
            continue
        t = f'{gene}'
        offset = 0
        text_labels.append(ax.annotate(t, (x[i], y[i]), xytext=(x[i] + offset, y[i] + offset), fontsize=10,
                                    ha='center', va='center', ))
        plt.scatter(x[i], y[i], s=8, zorder=3, facecolor = 'none', edgecolor='black')

    adjust_text(text_labels, ax=ax, lim=100,
                arrowprops= dict(arrowstyle='-', 
                            color='black', 
                            mutation_scale=6),
            )

    plt.xlabel("Count in Treg")
    plt.ylabel("Count in Tcon")
    plt.title("# Inter metadomains")
    plt.ylim([-10, 220])
    plt.xlim([-10, 220])
    return fig

def interval_to_str(interval):
    int1, int2 = interval.left, interval.right
    if np.isinf(int1):
        return f"≤{int(int2)}"
    if np.isinf(int2):
        return f"≥{int(int1)}"
    else:
        return f"{int(int1)}-{int(int2)}"

from adjustText import adjust_text
def inter_metadomains_vs_SE(all_inter_treg_metadomains, SE_count,):
    treg_inter = all_inter_treg_metadomains.sum(axis=1)
    n_ses = SE_count
    data = pd.DataFrame()
    data['# SEs'] = n_ses
    cuts = [0, 1, 5, 10, np.inf]
    data['# Inter metadomains'] = pd.cut(treg_inter, cuts).map(interval_to_str)
    fig, ax = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    sns.barplot(data=data, y='# SEs', x='# Inter metadomains', ax=ax, )
    return fig
    # scatter_with_pearson(x, y, s=8, ax=plt.gca(), zorder=3, color = 'skyblue', rasterized=True)

    # text_labels = []

    # genes_to_plot = set(list(np.where(y > 150)[0]) + list(np.where((x - y < -65))[0][::3]) + list(np.where(((x-y) > 20))[0]) + [9010])
    # for i in genes_to_plot:
    #     gene = get_name(i, ind_to_gene)
    #     if gene not in ['Tox', 'Ikzf2', 'Ikzf1', 'Pik3ip1', 'Lef1', 'Tgfbr3', 'Runx1', 'Lef1', 'Lta']:
    #         continue
    #     t = f'{gene}'
    #     offset = 0
    #     text_labels.append(ax.annotate(t, (x[i], y[i]), xytext=(x[i] + offset, y[i] + offset), fontsize=10,
    #                                 ha='center', va='center', ))
    #     plt.scatter(x[i], y[i], s=8, zorder=3, facecolor = 'none', edgecolor='black')

    # adjust_text(text_labels, ax=ax, lim=100,
    #             arrowprops= dict(arrowstyle='-', 
    #                         color='black', 
    #                         mutation_scale=6),
    #         )

    # plt.xlabel("Count in Treg")
    # plt.ylabel("Count in Tcon")
    # plt.title("# Inter metadomains")
    # plt.ylim([-10, 220])
    # plt.xlim([-10, 220])
    return fig




from adjustText import adjust_text
def inter_intra_metadomains(all_inter_treg_metadomains, all_intra_treg_metadomains, ind_to_gene, celltype):
    treg_inter = all_inter_treg_metadomains.sum(axis=1)
    treg_intra = all_intra_treg_metadomains.sum(axis=1)
    x, y = treg_inter, treg_intra

    fig, ax = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    scatter_with_pearson(x, y, plt.gca(), s=8, zorder=3, color = 'skyblue', rasterized=True)


    # print(np.where((x > 100) & (y > 40)))
    genes_to_plot = set(list(np.where(y > 40)[0]) + list(np.where(abs(x + y) > 40)[0]) + list(np.where(x > 50)[0]))
    genes_to_plot = [ 66, 3359, 6682, 278, 8848, 9010, 6121, 6155, 6154, 1885, 6682]
    for i in genes_to_plot:
        gene = get_name(i, ind_to_gene)
        # print(gene)
        if gene not in ['Ly96', 'Arpc1b', 'Adam17', 'Ikzf2', 'Runx1', 'Lta', 'Pik3ip1', 'Ikzf1',
                        'S100a11', 'Iah1']:
            continue
        t = f'{gene}'
        offset = 0
        # print(gene, x[i], y[i])
        text_labels.append(ax.annotate(t, (x[i], y[i]), xytext=(x[i] + offset, y[i] + offset), fontsize=10,
                                    ha='center', va='center', ))
        plt.scatter(x[i], y[i], s=8, zorder=3, facecolor = 'none', edgecolor='black')

    adjust_text(text_labels, ax=ax, lim=100,
                arrowprops = dict(arrowstyle='-', 
                            color='black', 
                            mutation_scale=6),
            )

    plt.xlabel(f"Inter, {celltype}")
    plt.ylabel(f"Intra, {celltype}")
    plt.title("# Intra vs # Inter metadomains")
    plt.ylim([-2.5, 100])
    plt.xlim([-5, 220])
    return fig


def inter_metadomain_barplot(all_inter_treg_metadomains):
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    plt.sca(axs)
    sns.histplot(all_inter_treg_metadomains.sum(axis=1), bins=100, rasterized=False,
                 ax=axs)
    plt.ylim([0, 500])
    plt.xlabel("# metadomains (Treg)")
    plt.gca().set_axisbelow(True)
    plt.title("# Inter metadomains")
    return fig


###

from plot_pvals import add_stat_annotation_boxplot_no_hue
def inter_vs_intra_h3k27ac_boxplot(all_intra_treg_metadomains, all_inter_treg_metadomains, 
                                   all_intra_metadomains, all_inter_metadomains, bw_val_df_all_250kb):
    all_mega = all_inter_metadomains.sum(axis=1) + all_intra_metadomains.sum(axis=1)
    delta = all_inter_treg_metadomains.sum(axis=1) - all_intra_treg_metadomains.sum(axis=1)
    
    good = all_mega > 10

    key = 'Treg H3K27ac'
    x, y = np.ravel(np.log2(1+bw_val_df_all_250kb[key])), all_intra_treg_metadomains.sum(axis=1)/all_mega
    x, y = np.ravel(np.log2(1+bw_val_df_all_250kb[key])), delta
    x, y = x[good], y[good]
    data = pd.DataFrame()
    data['x'] = x
    data['y'] = y
    data['xbins'] = pd.qcut(data['x'], np.linspace(0, 1, 11), labels=False)

    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)
    sns.boxplot(data=data, x='xbins', y='y', showfliers=False)
    
    # comparisons = list(zip(np.arange(10), np.arange(11)))
    comparisons = []
    comparisons.append([0, 4])
    comparisons.append([4, 9])
    add_stat_annotation_boxplot_no_hue(axs, data, 'xbins', 'y', np.arange(10), comparisons, 
                                       ymax=20, delta = 1, h = 2, log=False)    

    plt.xlabel("H3K27ac signal")
    plt.ylabel("# Inter - # Intra ")
    plt.title("Inter vs intra metadomains")
    return fig


import statsmodels.nonparametric.smoothers_lowess

def inter_vs_intra_h3k27ac_scatterplot(all_intra_treg_metadomains, all_inter_treg_metadomains, 
                                   all_intra_metadomains, all_inter_metadomains, bw_val_df_all_250kb):

    all_mega = all_inter_metadomains.sum(axis=1) + all_intra_metadomains.sum(axis=1)
    delta = all_inter_treg_metadomains.sum(axis=1) - all_intra_treg_metadomains.sum(axis=1)
    good = all_mega > 15

    key = 'Treg H3K27ac'
    x, y = np.ravel(np.log2(1+bw_val_df_all_250kb[key])), all_intra_treg_metadomains.sum(axis=1)/all_mega
    x, y = np.ravel(np.log2(1+bw_val_df_all_250kb[key])), delta

    # good[x < -2] = 0
    # y[y > 140] = 140
    x, y = x[good], y[good]
    smooth_x, smooth_y = statsmodels.nonparametric.smoothers_lowess.lowess(y, x, frac = .6).T

    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)
    plt.plot(smooth_x, smooth_y, color='red', zorder=3)
    plt.scatter(x, y, s = .1)
    plt.xlabel("Log2(H3K27ac)")
    plt.ylabel("#Inter - # Intra")
    plt.title(key)
    add_yaxis_labels('More Intra', 'More Inter', axs, fontsize=6)
    plt.ylim([-150, 150])
    plt.yticks([-75, 0, 75])
    return fig


def metadomain_density_barplot(metadomain_hub_freq_df):
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi=200)
    metadomain_hub_freq_df.plot.bar(stacked=True, color = ['salmon', 'lightgray', 'aqua', 'white'], ax=axs, linewidth=.1, edgecolor='black')
    axs.set_axisbelow(True)
    plt.xticks(rotation=0)
    plt.xlabel("Hub")
    plt.ylabel("Frac. of all \npossible metadomains")
    plt.title("Megaloop Density in Hub")
    leg = plt.legend(bbox_to_anchor=(0, 1), title='Megaloop in:', loc='upper left', frameon=False)
    leg.get_title().set_fontsize(6)
    plt.ylim([0, .50])
    return fig


def metadomain_between_hub_barplot(all_intra_metadomains, self, columns_to_names, row_colors, dpi = 200):
    metadomain_cluster_df = pd.DataFrame(all_intra_metadomains)
    metadomain_cluster_df['cluster'] = -1
    
    for u in [0, 4, 18]:
        inds = self.goodinds[self.merged_clustdict['all']==u]
        metadomain_cluster_df.loc[inds, 'cluster'] = u
    
    metadomain_cluster_df = metadomain_cluster_df.groupby('cluster').sum()
    metadomain_cluster_df = metadomain_cluster_df.T
    
    metadomain_cluster_df['cluster'] = -1
    
    for u in [0, 4, 18]:
        inds = self.goodinds[self.merged_clustdict['all']==u]
        metadomain_cluster_df.loc[inds, 'cluster'] = u
    
    metadomain_cluster_df = metadomain_cluster_df.groupby('cluster').sum()
    metadomain_cluster_df = (metadomain_cluster_df / metadomain_cluster_df.sum(axis=1)).T

    metadomain_cluster_df.columns = [columns_to_names.get(x, "All other bins") for x in metadomain_cluster_df.columns]
    metadomain_cluster_df.index = [columns_to_names.get(x, "All other bins") for x in metadomain_cluster_df.index]
    
    metadomain_cluster_df = metadomain_cluster_df.drop("All other bins")[['Active 1', 'Active 2', 'Repressive', 'All other bins']]
    
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi=dpi)
    metadomain_cluster_df.plot.bar(stacked=True, color = row_colors + ['lightgray'], ax=axs, linewidth=.1, edgecolor='black')
    axs.set_axisbelow(True)
    plt.xticks(rotation=0)
    plt.xlabel("Hub")
    plt.ylabel("# Metadomains")
    plt.title("Metadomains between Hubs")
    plt.legend(frameon=False, bbox_to_anchor=(1, 1))
    return fig

## Fig S15E
def treg_comp_violin_plot(my_treg_comp, self, columns_to_names):
    data = []
    for u in [0, 4, 18]:
        comp_v = pd.Series(my_treg_comp[self.goodinds[self.merged_clustdict['all']==u]])
        comp_v.index = [columns_to_names[u]]*len(comp_v)
        data.append(comp_v)
    comp_v = pd.Series(my_treg_comp)
    comp_v.index = ['All']*len(comp_v)
    data.append(comp_v)
    data = pd.concat(data, axis=0).reset_index()

    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 200)
    sns.violinplot(data=data, x='index',y=0,
                palette = ['lightgreen', 'green', 'orange'] + ['lightgray'])
    plt.xlabel("Hub")
    plt.ylabel("Compartment score")
    plt.yticks([0, -1, 1])

    add_stat_annotation_boxplot_no_hue(axs, data, 'index', 0, ['Active 1', 'Active 2', 'Repressive', 'All'],
                                    [['Active 2', 'Repressive'],
                                    ['Active 1', 'Active 2'], ], ymax = 1.35, h = .05
                                    )

    fig.savefig('./plots/paper/s15/E_compartment_violin.pdf', bbox_inches='tight')



def make_atac_enrichments_in_hub_quant(atac_peaks, self, all_ind_to_region, columns_to_names, row_colors_dict):
    vals = []
    all_atac_peak_df = atac_peaks.to_dataframe()
    all_lfcs = all_atac_peak_df['score']

    fig, axs = init_subplots_exact(3, 1, fgsz=(40*mm, 40*mm), dpi = 150,
                                  sharey=False)
    for c, u in enumerate(self.merged_inds_to_subset):
        plt.sca(axs[c])
        cluster_bedtools = add_chr_to_bedtool([all_ind_to_region[x] for x in self.goodinds[self.merged_clustdict['all']==u]])
        atac_peak_df = add_chr_to_bedtool(atac_peaks).intersect(cluster_bedtools, u=True).to_dataframe()
        lfcs = atac_peak_df['score']
        print(len(lfcs))
        sns.ecdfplot(lfcs, label=columns_to_names[u], 
                    color=row_colors_dict[columns_to_names[u]])
    
        sns.ecdfplot(all_lfcs, label='All peaks', 
                    color='black')
        add_pval_to_plot(lfcs, all_lfcs, plt.gca(),)
        
        plt.xlim([-2, 2])
        if c == 0:
            plt.ylabel("CDF")
        else:
            plt.ylabel("")
        plt.legend(bbox_to_anchor=(0, .95), loc='upper left', frameon=False)
        plt.xlabel("ATAC-seq LFC, peaks in hub")
        plt.title(columns_to_names[u])
    fig.savefig('./plots/paper/s15/quant_atac.pdf', bbox_inches='tight')    

def n_genes_per_hub(self, ind_to_gene, columns_to_names, row_colors_dict, my_treg_comp):
    n_genes = []
    baseline = [len(ind_to_gene.get(x, [])) for x in np.where(my_treg_comp > 0)[0]]

    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    for cluster in [0, 4, 18]:
        n_genes = [len(ind_to_gene.get(x, [])) for x in self.goodinds[self.merged_clustdict['all']==cluster]]
        p = format_pvalue(scipy.stats.ranksums(n_genes, baseline).pvalue)
        sns.ecdfplot(n_genes,
                    color = row_colors_dict[columns_to_names[cluster]],
                     linewidth=1,
                     label=f'{columns_to_names[cluster]}, p={p}',
                    )
    sns.ecdfplot(baseline, color = 'black', linewidth=1, label='All A comp.')
    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', frameon=False)
    plt.xlabel("# Genes per 250kb bin")    
    plt.title("# Genes in hub")
    plt.xlim([-1, 25])
    return fig

def gene_length_per_hub(self, gene_to_ind, ind_to_gene, columns_to_names, row_colors_dict, geneLengths):
    subgenes = geneLengths[geneLengths.index.isin(gene_to_ind.keys())]/1e3
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    for cluster in [0, 4, 18]:
        genes = np.unique(list(itertools.chain(*[ind_to_gene.get(x, []) for x in self.goodinds[self.merged_clustdict['all']==cluster]])))
        genes = pd.Series(genes)
        genes = genes[genes.isin(subgenes.index)]
        lengths = subgenes.loc[genes]
        p = format_pvalue(scipy.stats.ranksums(lengths, subgenes).pvalue)
        sns.ecdfplot(lengths, color = row_colors_dict[columns_to_names[cluster]],
                    linewidth = 1,
                    label = f'{columns_to_names[cluster]}, p={p}')
    sns.ecdfplot(subgenes, color = 'black', linewidth = 1, label = 'All')
    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', frameon=False)
    plt.xlabel("Gene Length (kb)")
    plt.title("Gene length in hub")
    plt.xlim([-1, 20])
    return fig

from scipy.stats import zscore
def plot_inter_metadomain_clustering_chip(bw_val_df_all_250kb, self, my_treg_comp):
    keysoi = [ 'Treg H3K4me3',  'Treg H3K4me1',  'Treg H3K27ac',  'Treg H3K27me3',  'Treg CTCF',  'Treg Smc1a', 
            ]
    # xticks = [
        # [0, 1, 2], 
    n = len(keysoi)
    v_df = pd.DataFrame()
    for c2, key in enumerate(keysoi):
        x = np.ravel(bw_val_df_all_250kb[key])
        c = np.asarray([-1]*len(x))
        for u in np.unique(self.clustdict['all']):
            idx = self.goodinds[self.clustdict['all']==u]
            c[idx] = u

        idx = ~np.isnan(x)
        x[idx] = zscore(x[idx])
        v_df[key] = x
        v_df['cluster'] = c

    idx = ~np.isnan(my_treg_comp)
    v_df['Compartment'] = my_treg_comp*np.nan
    v_df.loc[idx, 'Compartment'] = zscore(my_treg_comp[idx])

    # clustorder = list(self.cluster_to_subset_for_further_clustering)
    clustorder = list(np.unique(self.clustdict['all'][self.odict['all']]))

    v_df.loc[(v_df['cluster'] == -1), 'cluster'] = 'Others'
    data = v_df.dropna().melt(id_vars='cluster')
    data['cluster'] = data['cluster']
    data['variable'] = data['variable'].str.replace("Treg ", "")
    data = data[data['cluster']!=-1]
    fig, axs = init_subplots_exact(1, 1, dpi = 200, fgsz=(120*mm, 40*mm), yspace = 1.5, xspace=1.2, sharey=True)
    sns.boxplot(data=data,
                x='variable', y='value', hue='cluster',
                order = ['Compartment', 'H3K4me3', 
                        'H3K4me1', 'H3K27ac',
                        'Smc1a', 'CTCF', 
                'H3K27me3'
                ],
                hue_order = ['Others', ] + clustorder,
                palette = ['lightgray'] + [(sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20))[x] for x in clustorder],
                # orient='h', 
                showfliers=False, zorder=3,)
    plt.legend(bbox_to_anchor=(0, -.3), loc = 'center left', frameon=False, fontsize=6, ncol=13)
    plt.xlabel("")
    plt.ylabel("ChIP Enrichment")
    plt.title('ChIP Enrichment at Clusters')
    axs.set_axisbelow(True)
    return fig


def plot_inter_metadomain_clustering_chip2(bw_val_df_all_250kb, self, my_treg_comp):
    keysoi = [ 'Treg H3K4me3',  'Treg H3K4me1',  'Treg H3K27ac',  'Treg H3K27me3',  'Treg CTCF',  'Treg Smc1a', 
            ]
    # xticks = [
        # [0, 1, 2], 
    n = len(keysoi)
    v_df = pd.DataFrame()
    for c2, key in enumerate(keysoi):
        x = np.ravel(bw_val_df_all_250kb[key])
        c = np.asarray([-1]*len(x))
        for u in np.unique(self.clustdict['all']):
            idx = self.goodinds[self.clustdict['all']==u]
            c[idx] = u

        idx = ~np.isnan(x)
        x[idx] = zscore(x[idx])
        v_df[key] = x
        v_df['cluster'] = c

    idx = ~np.isnan(my_treg_comp)
    v_df['Compartment'] = my_treg_comp*np.nan
    v_df.loc[idx, 'Compartment'] = zscore(my_treg_comp[idx])

    # clustorder = list(self.cluster_to_subset_for_further_clustering)
    clustorder = list(np.unique(self.clustdict['all'][self.odict['all']]))

    v_df.loc[(v_df['cluster'] == -1), 'cluster'] = 'Others'
    data = v_df.dropna().melt(id_vars='cluster')
    data['cluster'] = data['cluster']
    data['variable'] = data['variable'].str.replace("Treg ", "")
    data = data[data['cluster']!=-1]
    fig, axs = init_subplots_exact(1, 1, dpi = 200, fgsz=(120*mm, 40*mm), yspace = 1.5, xspace=1.2, sharey=True)
    sns.boxplot(data=data,
                x='cluster', y='value', hue='variable',
                order = np.arange(0, 22),
                # hue_order = ['Others', ] + clustorder,
                # palette = ['lightgray'] + [(sns.color_palette('tab20b', n_colors=20) + sns.color_palette('tab20c', n_colors=20))[x] for x in clustorder],
                # orient='h', 
                showfliers=False, zorder=3,)
    plt.legend(bbox_to_anchor=(0, -.3), loc = 'center left', frameon=False, fontsize=6, ncol=13)
    plt.xlabel("")
    plt.ylabel("ChIP Enrichment")
    plt.title('ChIP Enrichment at Clusters')
    axs.set_axisbelow(True)
    return fig


def make_atac_enrichments_in_hub_with_pvals(self, all_ind_to_region, columns_to_names):
    atac_peaks = pbt.BedTool('./atac/processed/Treg_rest_vs_Tcon_rest_thresh=0.5.csv')
    vals = []
    for u in self.merged_inds_to_subset:
        cluster_bedtools = add_chr_to_bedtool([all_ind_to_region[x] for x in self.goodinds[self.merged_clustdict['all']==u]])
        atac_peak_df = add_chr_to_bedtool(atac_peaks).intersect(cluster_bedtools, u=True).to_dataframe()
        data_series = (np.sign(atac_peak_df['score']) * (atac_peak_df['strand'] < .05)).value_counts()
        data_series.name = columns_to_names[u]
        vals.append(data_series)

    all_atac_peak_df = atac_peaks.to_dataframe()
    data_series = (np.sign(atac_peak_df['score']) * (atac_peak_df['strand'] < .05)).value_counts()
    data_series.name = 'All'
    vals.append(data_series)

    data = pd.concat(vals, axis=1)
    data = data.loc[[0, -1, 1]]
    data.index = ['Common', 'Tcon', 'Treg']
    data = data.loc[['Treg', 'Common', 'Tcon']]

    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 200)
    data.div(data.sum(axis=0), axis=1).T.plot.bar(stacked=True, color=['red', 'lightgray', 'blue'], ax=axs)
    
    plt.gca().set_axisbelow(True)
    plt.legend(bbox_to_anchor=(1,1))
    plt.xticks(rotation=0)
    plt.ylabel("Fraction")
    plt.title("Fraction differential ATAC peaks")
    
    for i, column in enumerate(data.columns):
        for j, row in enumerate(data.index):
            # Perform Fisher's Exact Test
            table = np.array([[data.loc[row, column], data[column].sum() - data.loc[row, column]],
                              [data.loc[row, "All"], data['All'].sum() - data.loc[row, "All"]]])
            # print(column, row, table)
            _, p_value = scipy.stats.fisher_exact(table)
    
            # Outline the box if significant
            if p_value < 0.05:
                print(column, row)
                rect = axs.patches[j*4 + i]
                rect.set_edgecolor('black')
                rect.set_linewidth(.5)
    fig.savefig('./plots/paper/s15/atac_binary.pdf', bbox_inches='tight')
    

def housekeeping_barplots(self, gene_to_ind, ind_to_gene, row_colors):
    all_genes = set(list(gene_to_ind.keys()))
    housekeeping_genes2 = set(list(pd.read_csv('./housekeeping_genes/Bulk Cluster Ubiquitous.txt', sep='\t', header=None).T.dropna()[0].values))
    housekeeping_genes2 = housekeeping_genes2 & all_genes
    
    housekeeping_genes = set(list(pd.read_csv('./housekeeping_genes/Housekeeping_GenesMouse.csv', sep=';')['Gene']))
    housekeeping_genes = housekeeping_genes & all_genes
    
    basep = len(all_genes & housekeeping_genes) / len(all_genes)
    
    fig, axs = init_subplots_exact(1, 1, fgsz=(60*mm, 40*mm), dpi = 150)
    colors = list(sns.color_palette('tab20b')) + list(sns.color_palette('tab20c'))
    ps = []
    clusters = np.unique(self.clustdict['all'])
    
    sig = []
    for u in clusters:
        genes = set(list(itertools.chain(*[ind_to_gene.get(x, []) for x in self.goodinds[self.clustdict['all']==u]])))
        p = len(genes & housekeeping_genes)/len(genes)
        stat, pval = scipy.stats.fisher_exact([[len(genes & housekeeping_genes), len(genes)],
                                        [len(all_genes & housekeeping_genes), len(all_genes), ]]
                                       )
        if pval < .05:
            pval = format_pvalue(pval)
            plt.text(u, p*1.02, pval, fontsize=6, ha='center', rotation=90, va='bottom')
        sig.append(pval)
        ps.append(p)
    
    plt.bar(clusters, ps, color = colors[:len(ps)])
    plt.ylim([0, .35])
    plt.ylabel("Fraction")
    plt.xlabel("Cluster (unmerged)")
    plt.axhline(basep, color='black', linestyle='--')
    plt.title("Fraction housekeeping genes")
    fig.savefig("./plots/paper/s15/n_housekeeping_cluster.pdf", bbox_inches='tight')
    
    basep = len(all_genes & housekeeping_genes) / len(all_genes)
    
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    colors = list(sns.color_palette('tab20b')) + list(sns.color_palette('tab20c'))
    ps = []
    clusters = np.unique([0, 4, 18])
    
    sig = []
    for c, u in enumerate(clusters):
        genes = set(list(itertools.chain(*[ind_to_gene.get(x, []) for x in self.goodinds[self.merged_clustdict['all']==u]])))
        p = len(genes & housekeeping_genes)/len(genes)
        stat, pval = scipy.stats.fisher_exact([[len(genes & housekeeping_genes), len(genes)],
                                        [len(all_genes & housekeeping_genes), len(all_genes), ]]
                                       )
    
    
        pval = format_pvalue(pval)
        plt.text(c, p*1.02, f'{pval}\nn={len(genes & housekeeping_genes)}/{len(genes)}', fontsize=6, ha='center', va='bottom')
        sig.append(pval)
        ps.append(p)
    xs = np.arange(3)
    plt.bar(xs, ps, color=row_colors)
    plt.ylim([0, .35])
    plt.ylabel("Fraction")
    plt.xlabel("Cluster (unmerged)")
    plt.axhline(basep, color='black', linestyle='--')
    plt.title("Fraction housekeeping genes")
    fig.savefig("./plots/paper/s15/n_housekeeping_hub.pdf", bbox_inches='tight')    

### Figure S20
import matplotlib.patches as mpatches
from collections import Counter
def drop_duplicate_motif_families(families, k=3):
    to_keep = []
    counter = Counter()
    for factor, family in families.items():
        counter[family] += 1
        if counter[family] < k:
            to_keep.append(factor)
    return to_keep

def make_motif_enrichment_plot(self, _250kb_hub_annotations, all_ind_to_region, atac_peak_motif_counts, motif_metadata, columns_to_names):    
    motif_result_df_dict = {}
    for u in ['matched_A compartment', 'Constitutive', 'Dynamic', 'Repressive']:
        cluster_bedtools = add_chr_to_bedtool([all_ind_to_region[x] for x in _250kb_hub_annotations.index[_250kb_hub_annotations['Hub']==u]
                                               ])
        peak_idx = bedtool_to_index(index_to_bedtool(atac_peak_motif_counts.index).intersect(cluster_bedtools, u=True))
        vals = []
        xs = atac_peak_motif_counts.loc[peak_idx]
        ys = atac_peak_motif_counts.drop(peak_idx)
        for col in atac_peak_motif_counts:
            x = xs[col]
            y = ys[col]
            _, p = scipy.stats.ranksums(x, y)
            delta = np.log2(x.mean() / y.mean())
            vals.append([col, p, delta])
        motif_result_df = pd.DataFrame(vals, columns = ['factor', 'p', 'delta'])
        motif_result_df_dict[u] = motif_result_df

    keys_to_add = []
    motif_data = pd.DataFrame()
    for key in motif_result_df_dict:
        data = motif_result_df_dict[key]
        up_inds = data[(data['p'] < 1e-4) & (data['delta'] > 0)].sort_values('p')['factor'].iloc[:40]
        down_inds = data[(data['p'] < 1e-4) & (data['delta'] < 0)].sort_values('p')['factor'].iloc[:40]
        up_inds = drop_duplicate_motif_families(motif_metadata.loc[up_inds, 'coarse_parsed_family'], k=3)
        down_inds = drop_duplicate_motif_families(motif_metadata.loc[down_inds, 'coarse_parsed_family'], k=3)
        up_inds = up_inds[:10]
        down_inds = down_inds[:10]

        keys_to_add.extend(list(up_inds))
        keys_to_add.extend(list(down_inds))

    keys_to_add = np.unique(keys_to_add)

    for key in motif_result_df_dict:
        data = motif_result_df_dict[key].set_index('factor')
        motif_data[key] = (data['delta'] * (data['p'] < .2)).loc[keys_to_add]


    coarse_parsed_families = motif_metadata.loc[keys_to_add, 'coarse_parsed_family'].unique()
    colors = sns.color_palette('hls', n_colors=len(coarse_parsed_families))
    family_color_dict = dict(zip(coarse_parsed_families, colors))


    motif_data.index.name = 'Factor'

    row_colors = motif_metadata.loc[keys_to_add, 'coarse_parsed_family'].apply(family_color_dict.get)
    row_colors.name = 'Motif Family'

    g = sns.clustermap(motif_data, vmin=-1, vmax=1, cmap='coolwarm',
                row_colors=row_colors, figsize=(4, 6), yticklabels=True,
                linewidth=.05,)

    # Function to create legend
    def add_legend(ax, labels, colors):
        patches = [mpatches.Patch(color=color, label=label) for label, color in zip(labels, colors)]
        legend = ax.legend(handles=patches, loc='upper left', bbox_to_anchor=(1.2, 1), ncol=2, title='Motif Family')
        return legend

    # Extract unique labels and their corresponding colors
    unique_labels = family_color_dict.keys()
    unique_colors = [family_color_dict[label] for label in unique_labels]

    # Add legend to the plot
    add_legend(g.ax_heatmap, unique_labels, unique_colors)
    g.savefig('./plots/paper/s20/motifs_at_hubs.pdf', bbox_inches='tight')

    # fig.savefig('./plots/motif_enrichment/metadomains.pdf', bbox_inches='tight')
    return motif_result_df_dict


def prepare_cusanovich_data():
    cell_metadata = pd.read_csv('./mouse_atac_atlas/cell_metadata.txt', sep='\t')
    cell_metadata['fine_cluster'] = cell_metadata['cluster'].astype(str) + "_" + cell_metadata['subset_cluster'].astype(str)
    matrix = scipy.io.mmread('./mouse_atac_atlas/matrix.mtx.gz').tocsr()

    pseudobulk_data = pd.DataFrame()
    factor = 'cluster'
    for u in cell_metadata[factor].unique():
        cell_idx = cell_metadata[factor]==u
        pseudobulk_data[u] = np.ravel(matrix[:, cell_idx].mean(axis=1))
        print("Done with", u)

    pseudobulk_data.columns = cell_metadata[factor].unique()
    pseudobulk_data.to_csv('./mouse_atac_atlas/pseudobulk_data.csv')
    return pseudobulk_data

def cusanovich_atac_clustering(self, _250kb_hub_annotations):

    def zscore_normalize(data, axis=0):
        """
        Normalize the DataFrame using Z-score normalization along the specified axis.

        Parameters:
        data (pd.DataFrame): The input DataFrame.
        axis (int): The axis along which to compute the Z-scores. 
                    0 for columns and 1 for rows.

        Returns:
        pd.DataFrame: The normalized DataFrame.
        """
        if axis == 0:
            mean = np.mean(data.values, axis=0)
            std = np.std(data.values, axis=0)
            zscore_data = (data.values - mean) / std
        elif axis == 1:
            mean = np.mean(data.values, axis=1).reshape(-1, 1)
            std = np.std(data.values, axis=1).reshape(-1, 1)
            zscore_data = (data.values - mean) / std
        else:
            raise ValueError("Axis must be 0 (columns) or 1 (rows).")
        
        return pd.DataFrame(zscore_data, index=data.index, columns=data.columns)
    
    factor = 'cluster'
    pseudobulk_data = pd.read_csv('./mouse_atac_atlas/pseudobulk_data.csv', sep='\t', index_col=0)
    cell_metadata = pd.read_csv('./mouse_atac_atlas/cell_metadata.txt', sep='\t')
    cell_metadata['fine_cluster'] = cell_metadata['cluster'].astype(str) + "_" + cell_metadata['subset_cluster'].astype(str)

    zscore_pseudobulk_data_0 = zscore_normalize(pseudobulk_data, axis=0)
    # zscore_pseudobulk_data_1 = zscore_normalize(pseudobulk_data, axis=1)

    peak2loc = {}
    for u in ['matched_A compartment', 'Constitutive', 'Dynamic', 'Repressive']:
        cluster_bedtools = add_chr_to_bedtool([self.all_ind_to_region[x] for x in 
                                               _250kb_hub_annotations.index[_250kb_hub_annotations['Hub']==u]
                                               ])

        peak_idx = np.where(index_to_bedtool(np.ravel(pd.read_csv('./mouse_atac_atlas/atac_matrix.binary.qc_filtered.peaks.txt', header=None).values)
                        ).intersect(cluster_bedtools, c=True).to_dataframe()['name']>0)[0]
        peak2loc[u] = peak_idx


    data = []
    us = ['matched_A compartment', 'Constitutive', 'Dynamic', 'Repressive']
    for u in us:
        peak_idx = peak2loc[u]
        submat = np.ravel(zscore_pseudobulk_data_0.loc[peak_idx].mean(axis=0))
        submat = pd.DataFrame(submat)
        data.append(submat)

    data = pd.concat(data, axis=1)
    data.columns = us
    data.index = np.asarray(zscore_pseudobulk_data_0.columns).astype(int)

    cellfreqs = cell_metadata[[factor, 'tissue']].value_counts().unstack().fillna(0)
    cellfreqs = (cellfreqs / cellfreqs.sum(axis=0))
    cluster2cellfreq = pd.Series(cellfreqs.columns[np.argmax(cellfreqs, axis=1)], index=cellfreqs.index).loc[data.index]
    data.index = cluster2cellfreq
    g = sns.clustermap(data, cmap='coolwarm',
                figsize=(3, 4), yticklabels=True,
                col_colors = ['lightgray', 'lightgreen', 'green', 'orange', ],)
    g.savefig('./plots/paper/fig3/cusanovich_atac.pdf', bbox_inches='tight')
    return g.fig






### Figure S21
from tad_functions import make_order2_and_cluster
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

def metadomain_vs_oe_clustering(self, sep_oe_mat_treg, sep_oe_mat_tcon, bw_val_df_all_250kb, mm=1/2.54/10):
    np.random.seed(1)
    randinds = self.goodinds
    # randinds = np.random.choice(np.where(~np.isnan(my_treg_comp))[0], 
                                # size=len(self.goodinds),
                                # replace=False,
                            # )
    mtmp = (sep_oe_mat_tcon[randinds, :][:, randinds].copy() + sep_oe_mat_treg[randinds, :][:, randinds].copy())/2
    mtmp[np.isnan(mtmp)] = 0
    mtmp[np.isinf(mtmp)] = -10
    pca = PCA(n_components=20)
    pca = pca.fit(mtmp)
    components = pca.transform(mtmp)
    tmp_o, tmp_clusts = make_order2_and_cluster(components, n_clusters=24, method='ward', metric='euclidean')

    # def filter_low_clust(t, co = 10):
    #     z = t.value_counts('clust')
    #     return t[t['clust'].isin(z.index[z > co])]
        
    df = pd.DataFrame()
    for key in ['Treg H3K27ac',  'Treg H3K4me3', 'Treg H3K4me1']:
        df[key] = np.ravel(self.bw_val_dict[key])

    df['SE'] = self.SE_count_dict['common_SE'] > 0
    treal = df.copy()
    treal['clust'] = self.clustdict['all']

    t2 = pd.DataFrame()
    for key in ['Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1']:
        # t2[key] = np.ravel(self.bw_val_dict[key])
        t2[key] = np.ravel(bw_val_df_all_250kb[key][randinds].values)
    t2['clust'] = tmp_clusts
    t2['SE'] = self.SE_count_dict['common_SE']>0
    t2['ind'] = randinds
    cols = treal.columns[:-1]

    # for c, key in enumerate(['SE', 'Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1']):
        # labels = ['Clustering of metadomains', 'Clustering of O/E']
        # for df, label in zip([treal, t2], labels):
            # score = silhouette_score(df[key].values.reshape(-1, 1), np.random.permutation(df['clust']))
            # print(key, label, scipy.stats.f_oneway(df[key], df['clust']))
            # print(key, label)
            # print("__")

    n = len(cols)
    fig, axs = init_subplots_exact(n, 2, fgsz=(40*mm, 40*mm), dpi = 100, sharex=True, yspace=1.4, xspace=1.4)
    for c, key in enumerate(['SE', 'Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1']):
        plt.sca(axs[c])
        labels = ['Clustering of metadomains', 'Clustering of O/E']
        colors = ['purple', 'skyblue']

        for df, label, color in zip([treal, t2], labels, colors):
            df = df
            if 'SE' not in key:
                df[key] = scipy.stats.zscore(1+df[key])
                plt.ylabel("Zscore, ChIP Signal")
            else:
                plt.ylabel("Count")
            reorder_dict = dict(zip(df.groupby('clust').mean()[key].sort_values(ascending=False).reset_index()['clust'],
                                    df.groupby('clust').mean()[key].sort_values(ascending=False).reset_index().index)
                            )
            df['new_clust']  = df['clust'].apply(reorder_dict.get)
            sns.lineplot(data=df, x = 'new_clust', y = key, label = label,  c=color, 
                        marker='o', markersize=0,
                        )
            sns.scatterplot(data=df.groupby('new_clust').mean()[key].reset_index(), x = 'new_clust', y = key, label = None,  c=color, 
                        s=df['new_clust'].value_counts().sort_index()/3,
                        )
                
        for i in range(2):
            real, just_hic = treal[key][treal['new_clust']==i], t2[key][t2['new_clust']==i]
            if key == 'SE':
                p = scipy.stats.fisher_exact([[(real>0).sum(), (real==0).sum()], [(just_hic>0).sum(), (just_hic==0).sum()]])[1]                
            else:
                p = scipy.stats.ranksums(real, just_hic)[1]
            if (p < .05) and (key == 'SE') and (real.mean() > just_hic.mean()):
                plt.text(i, .8, '*', color='purple', ha='center')
            elif (p < .05) and (key != 'SE') and (real.mean() > just_hic.mean()):
                plt.text(i, 2.5, '*', color='purple', ha='center')
            elif (p < .05) and (key == 'SE') and (real.mean() < just_hic.mean()):
                plt.text(i, .8, '*', color='lightblue', ha='center')
            elif (p < .05) and (key != 'SE') and (real.mean() < just_hic.mean()):
                plt.text(i, 2.5, '*', color='lightblue', ha='center')
        

        plt.title(key)
        plt.legend(frameon=False, handlelength=1)
        plt.xlabel("Rank")
    fig.savefig('./plots/paper/s21/hic_oe_clustering.pdf', bbox_inches='tight')
    return t2


def make_df_from_inds(SE_count_dict, bw_val_df_all_250kb, inds, all_ind_to_region, cols, my_treg_comp):
    df = pd.DataFrame()
    for key in cols:
        df[key] = np.ravel(bw_val_df_all_250kb[key])[inds]
    df['SE'] = (SE_count_dict['common_SE'] > 0)[inds]
    df['ind'] = inds
    df['chrom'] = [all_ind_to_region[x][0] for x in inds]
    df['comp'] = my_treg_comp[inds]
    return df


def metadomain_vs_random_bin_clustering(self, sep_oe_mat_treg, sep_oe_mat_tcon, bw_val_df_all_250kb, my_treg_comp,
                                        SE_count_dict, all_ind_to_region, mm=1/2.54/10, method='random',
                                        ):
    np.random.seed(0)
    if method == 'random':
        randinds = np.random.choice(np.where(~np.isnan(my_treg_comp))[0], 
                                size=len(self.goodinds),
                                replace=False,
                            )
    elif method == 'same':
        randinds = self.goodinds
    elif method == 'acomp':
        randinds = np.argsort(-my_treg_comp)[:len(self.goodinds)]
        print("Using acomps")
    elif method == 'acomp_no_overlap':
        randinds = np.argsort(-my_treg_comp)
        randinds = randinds[~np.isin(randinds, self.goodinds)][:len(self.goodinds)]
    mtmp = (sep_oe_mat_tcon[randinds, :][:, randinds].copy() + sep_oe_mat_treg[randinds, :][:, randinds].copy())/2
    mtmp[np.isnan(mtmp)] = 0
    mtmp[np.isinf(mtmp)] = -10
    pca = PCA(n_components=20)
    pca = pca.fit(mtmp)
    components = pca.transform(mtmp)
    tmp_o, tmp_clusts = make_order2_and_cluster(components, n_clusters=24, method='ward', metric='euclidean')

    # def filter_low_clust(t, co = 10):
    #     z = t.value_counts('clust')
    #     return t[t['clust'].isin(z.index[z > co])]
    
    cols = ['Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1', 'Treg H3K27me3', 'Treg ATAC', 
            #'Treg CTCF', 'Treg Stat5', 'Treg Satb1'
            ]
    treal = make_df_from_inds(SE_count_dict, bw_val_df_all_250kb, self.goodinds, all_ind_to_region, cols,
                                    my_treg_comp)
    treal['clust'] = self.clustdict['all']

    trand = make_df_from_inds(SE_count_dict, bw_val_df_all_250kb, randinds, all_ind_to_region, cols,
                                    my_treg_comp)
    trand['clust'] = tmp_clusts
    
        
    # for c, key in enumerate(['SE', 'Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1', 'Treg H3K27me3']):
        # labels = ['Clustering of metadomains', 'Clustering of O/E']
        # for df, label in zip([treal, t2], labels):
            # print(key, label, scipy.stats.f_oneway(df[key], df['clust']))
            # print("__")
    cols = ['SE'] + cols
    n = len(cols)
    fig, axs = init_subplots_exact(n, 2, fgsz=(40*mm, 40*mm), dpi = 100, sharex=True, yspace=1.4, xspace=1.4)
    for c, key in enumerate(cols):
        plt.sca(axs[c])
        labels = ['Clustering of metadomains', 'Clustering of O/E']
        colors = ['purple', 'skyblue']

        for df, label, color in zip([treal, trand], labels, colors):
            df = df
            if 'SE' not in key:
                df[key] = df[key]
                plt.ylabel("Normalized ChIP Signal")
            else:
                plt.ylabel("Count")
            reorder_dict = dict(zip(df.groupby('clust')[key].mean().sort_values(ascending=False).reset_index()['clust'],
                                    df.groupby('clust')[key].mean().sort_values(ascending=False).reset_index().index)
                            )
            df['new_clust']  = df['clust'].apply(reorder_dict.get)
            sns.lineplot(data=df, x = 'new_clust', y = key, label = label,  c=color, 
                        marker='o', markersize=0,
                        )
            sns.scatterplot(data=df.groupby('new_clust')[key].mean().reset_index(), x = 'new_clust', y = key, label = None,  c=color, 
                        s=df['new_clust'].value_counts().sort_index()/3,
                        )
                
        for i in range(2):
            real, just_hic = treal[key][treal['new_clust']==i], trand[key][trand['new_clust']==i]
            if key == 'SE':
                p = scipy.stats.fisher_exact([[(real>0).sum(), (real==0).sum()], [(just_hic>0).sum(), (just_hic==0).sum()]])[1]                
            else:
                p = scipy.stats.ranksums(real, just_hic)[1]
            
            y = plt.ylim()[1] - .1
            if (p < .05) and (key == 'SE') and (real.mean() > just_hic.mean()):
                plt.text(i, .8, '*', color='purple', ha='center')
            elif (p < .05) and (key != 'SE') and (real.mean() > just_hic.mean()):
                plt.text(i, y, '*', color='purple', ha='center')
            elif (p < .05) and (key == 'SE') and (real.mean() < just_hic.mean()):
                plt.text(i, .8, '*', color='lightblue', ha='center')
            elif (p < .05) and (key != 'SE') and (real.mean() < just_hic.mean()):
                plt.text(i, y, '*', color='lightblue', ha='center')
        

        plt.title(key)
        plt.legend(frameon=False, handlelength=1)
        plt.xlabel("Rank")
    fig.savefig(f'./plots/paper/s21/hic_oe_clustering_method={method}.pdf', bbox_inches='tight')
    return treal, trand

def metadomain_vs_random_clustering(self, sep_oe_mat_treg, sep_oe_mat_tcon, bw_val_df_all_250kb, mm=1/2.54/10):
    np.random.seed(1)
    randinds = self.goodinds
    # randinds = np.random.choice(np.where(~np.isnan(my_treg_comp))[0], 
                                # size=len(self.goodinds),
                                # replace=False,
                            # )
    mtmp = (sep_oe_mat_tcon[randinds, :][:, randinds].copy() + sep_oe_mat_treg[randinds, :][:, randinds].copy())/2
    mtmp = np.random.permutation(mtmp)
    mtmp = np.random.permutation(mtmp.T)
    mtmp[np.isnan(mtmp)] = 0
    mtmp[np.isinf(mtmp)] = -10
    pca = PCA(n_components=20)
    pca = pca.fit(mtmp)
    components = pca.transform(mtmp)
    tmp_o, tmp_clusts = make_order2_and_cluster(components, n_clusters=24, method='ward', metric='euclidean')

    # def filter_low_clust(t, co = 10):
    #     z = t.value_counts('clust')
    #     return t[t['clust'].isin(z.index[z > co])]
        
    df = pd.DataFrame()
    for key in ['Treg H3K27ac',  'Treg H3K4me3', 'Treg H3K4me1']:
        df[key] = np.ravel(self.bw_val_dict[key])

    df['SE'] = self.SE_count_dict['common_SE'] > 0
    treal = df.copy()
    treal['clust'] = self.clustdict['all']

    t2 = pd.DataFrame()
    for key in ['Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1']:
        # t2[key] = np.ravel(self.bw_val_dict[key])
        t2[key] = np.ravel(bw_val_df_all_250kb[key][randinds].values)
    t2['clust'] = tmp_clusts
    t2['SE'] = self.SE_count_dict['common_SE']>0


    cols = treal.columns[:-1]
    n = len(cols)
    fig, axs = init_subplots_exact(n, 2, fgsz=(40*mm, 40*mm), dpi = 100, sharex=True, yspace=1.4, xspace=1.4)
    for c, key in enumerate(['SE', 'Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1']):
        plt.sca(axs[c])
        labels = ['Clustering of metadomains', 'Clustering of O/E']
        colors = ['purple', 'skyblue']

        for df, label, color in zip([treal, t2], labels, colors):
            df = df
            if 'SE' not in key:
                df[key] = scipy.stats.zscore(1+df[key])
                plt.ylabel("Zscore, ChIP Signal")
            else:
                plt.ylabel("Count")
            reorder_dict = dict(zip(df.groupby('clust').mean()[key].sort_values(ascending=False).reset_index()['clust'],
                                    df.groupby('clust').mean()[key].sort_values(ascending=False).reset_index().index)
                            )
            df['new_clust']  = df['clust'].apply(reorder_dict.get)
            sns.lineplot(data=df, x = 'new_clust', y = key, label = label,  c=color, 
                        marker='o', markersize=0,
                        )
            sns.scatterplot(data=df.groupby('new_clust').mean()[key].reset_index(), x = 'new_clust', y = key, label = None,  c=color, 
                        s=df['new_clust'].value_counts().sort_index()/3,
                        )
                
        for i in range(2):
            real, just_hic = treal[key][treal['new_clust']==i], t2[key][t2['new_clust']==i]
            if key == 'SE':
                p = scipy.stats.fisher_exact([[(real>0).sum(), (real==0).sum()], [(just_hic>0).sum(), (just_hic==0).sum()]])[1]
                print(p)
                
            else:
                p = scipy.stats.ranksums(real, just_hic)[1]
            if (p < .05) and (key == 'SE') and (real.mean() > just_hic.mean()):
                plt.text(i, .8, '*', color='purple', ha='center')
            elif (p < .05) and (key != 'SE') and (real.mean() > just_hic.mean()):
                plt.text(i, 2.5, '*', color='purple', ha='center')
            elif (p < .05) and (key == 'SE') and (real.mean() < just_hic.mean()):
                plt.text(i, .8, '*', color='lightblue', ha='center')
            elif (p < .05) and (key != 'SE') and (real.mean() < just_hic.mean()):
                plt.text(i, 2.5, '*', color='lightblue', ha='center')
        

        plt.title(key)
        plt.legend(frameon=False, handlelength=1)
        plt.xlabel("Rank")
    fig.savefig('./plots/paper/s21/hic_random_clustering.pdf', bbox_inches='tight')



   

import itertools
def baseline_scrna_correlation(gene_dict, ind_to_gene, self, tregs_rna_corr, row_colors_dict, columns_to_names):
    fig, axs = init_subplots_exact(3, 1, dpi=200, fgsz=(40*mm, 40*mm),
                                sharey=True)
    for c, col in enumerate([0, 4, 18]):
        genes_in_cluster = list(itertools.chain(*[ind_to_gene.get(x, []) for x in self.goodinds[self.merged_clustdict['all']==col]]))
        genes_in_cluster = pd.Series(genes_in_cluster).str.lower().values
        genes_in_cluster = genes_in_cluster[np.isin(genes_in_cluster, tregs_rna_corr.index)]
        genes_not_in_cluster = tregs_rna_corr.index.drop(genes_in_cluster)
        plt.sca(axs[c])

        corrs_treg = np.ravel(np.triu(tregs_rna_corr.loc[genes_in_cluster, genes_in_cluster].values, k = 1))
        corrs_treg = corrs_treg[corrs_treg != 0]
        sns.ecdfplot(corrs_treg, label=f'Hub <-> Hub; n={len(corrs_treg)}', color=row_colors_dict[columns_to_names[col]])

        corrs_treg_nonhub = np.ravel(np.triu(tregs_rna_corr.loc[genes_in_cluster, genes_not_in_cluster].values, k = 1))
        corrs_treg_nonhub = corrs_treg_nonhub[corrs_treg_nonhub != 0]
        sns.ecdfplot(corrs_treg_nonhub, label=f'Hub <-> Non-hub; n={len(corrs_treg_nonhub)}', linestyle='--', color=row_colors_dict[columns_to_names[col]])
        print(nonan_test(corrs_treg, corrs_treg_nonhub))

        plt.legend(frameon=False)
        plt.title(columns_to_names[col])
        if c > 0:
            axs[c].set_ylabel("")
        plt.xlabel("Pearson R")
        plt.xlim([-.3, .3])    
    return fig

def differential_scrna_correlation(gene_dict, ind_to_gene, self, tregs_rna_corr, treg_precursor_rna_corr, 
hub_pileup_pval_df_250kb, hub_pileup_stat_df_250kb, columns_to_names,
                                   diff_hic_pco, diff_hic_statco,
                                   ):
    treg_sets = []
    tcon_sets = []
    geneset = set(gene_dict['Resting'].index)
    fig, axs = init_subplots_exact(4, 2, dpi=200, fgsz=(40*mm, 40*mm),
                                sharey=True)
    for c, col in enumerate([0, 4]):
        treg_up = hub_pileup_pval_df_250kb.index[(hub_pileup_pval_df_250kb[col] < diff_hic_pco) & (hub_pileup_stat_df_250kb[col] > diff_hic_statco)]
        ns_bins = hub_pileup_pval_df_250kb.index[((hub_pileup_pval_df_250kb[col] > .5)) & (hub_pileup_stat_df_250kb[col].abs() < .05)]
        
        treg_sets.append(set(treg_up))

        name = columns_to_names[col].replace('\n', " ")
        set_dict = {}
        set_dict[f'Treg-up Hi-C'] = set(treg_up)
        set_dict['NS'] = set(ns_bins)
        
        kwarg_dict = {}
        kwarg_dict[f'Treg-up Hi-C'] = {}
        kwarg_dict[f'NS'] = {'linestyle' : '--'}
        kwarg_dict[f'Treg'] = {'color': 'red'}
        kwarg_dict[f'Prec.'] = {'color': 'gray'}

        
        genes_in_cluster = list(itertools.chain(*[ind_to_gene.get(x, []) for x in self.goodinds[self.merged_clustdict['all']==col]]))
        genes_in_cluster = pd.Series(genes_in_cluster).str.lower().values
        genes_in_cluster = genes_in_cluster[np.isin(genes_in_cluster, tregs_rna_corr.index)]
        

        for c2, (celltype, corr) in enumerate({'Treg' : tregs_rna_corr, 'Prec.' : treg_precursor_rna_corr}.items()):
            treg_val_dict = {}
            for name, indsoi in set_dict.items():
                genes = []
                for x in indsoi:
                    genes.extend(ind_to_gene.get(x, []))
                genes = arr(list(set(genes).intersection(geneset)))
                genes = pd.Series(genes).str.lower().values
                genes = genes[np.isin(genes, tregs_rna_corr.index)]
        
                genes_in_just_cluster = genes_in_cluster[~np.isin(genes_in_cluster, genes)]
                corrs_v = corr.loc[genes, genes_in_just_cluster]
                plt.sca(axs[c+c2*2])
                sns.ecdfplot(corrs_v.mean(axis=1), label=f'{celltype} RNA; {name} bins', **kwarg_dict[name], **kwarg_dict[celltype])
                
                treg_val_dict[f'{name}'] = corrs_v.mean(axis=1)

            p = format_pvalue(scipy.stats.ranksums(treg_val_dict['Treg-up Hi-C'], treg_val_dict['NS'])[1])
            plt.text(1, 0, p, transform=plt.gca().transAxes, va='bottom', ha='right')
        plt.sca(axs[c])
        plt.title(columns_to_names[col])

        for a in axs:
            plt.sca(a)
            plt.legend(frameon=False, handlelength=1, loc='upper left')
        if c > 0:
            axs[c].set_ylabel("")
        plt.xlabel("Pearson R")
    return fig
    

### Figure S24
def superenhancer_metadomain_score_cdf_plot(hub_pileup_stat_df_250kb, SE_treg_count, SE_tcon_count, SE_common_count, columns_to_names):
    fig, axs = init_subplots_exact(3, 1, fgsz=(40*mm, 40*mm), dpi = 100, xspace=1.4)
    for c, col in enumerate([0, 4, 18]):
        plt.sca(axs[c])
        sns.ecdfplot(hub_pileup_stat_df_250kb.loc[np.where(SE_treg_count>0)[0]][col], color='red', label='Treg SE')
        sns.ecdfplot(hub_pileup_stat_df_250kb.loc[np.where(SE_tcon_count>0)[0]][col], color='blue', label='Tcon SE')
        sns.ecdfplot(hub_pileup_stat_df_250kb.loc[np.where(SE_common_count>0)[0]][col], color='lightgray', label='Shared SE')
        add_xaxis_labels("Tcon", "Treg", plt.gca(), fontsize=8)
        plt.legend()
        plt.xlabel("Megaloop score (MS)")
        plt.title(columns_to_names[col])
    return fig




### Figure S28

def metaloop_h3k27ac_boxplot(all_h3k27ac, metaloop_anchors):
    l2 = add_chr_to_bedtool(pbt.BedTool.from_dataframe(all_h3k27ac)).intersect(metaloop_anchors, c=True).to_dataframe(header=None)

    # l2['quantiles'] = pd.cut(l2['thickStart'], [-np.inf, .5, 1.5, 2.5, np.inf], labels=False)
    l2['quantiles'] = l2['thickStart'].clip(0, 4)
    xs = []
    ys = [] 
    ns = []
    colors = sns.color_palette("coolwarm", n_colors = 5)
    # for c, (clow, chigh) in enumerate(cutoffs):
    # names = l[3][(l['score'] >= clow) & (l['score'] <= chigh)]
    x = l2['quantiles']
    y = np.log2(l2['name'])
    data = pd.concat([x, y], axis=1)
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)
    sns.boxplot(data=data, x='quantiles', y='name', showfliers=False)
    plt.ylabel("Log2(H3K27ac)\n(peaks in metaloop anchor)")
    plt.xlabel("# Metaloops")

    labels = [int(x.get_text()) for x in plt.gca().get_xticklabels()]
    labels[-1] = '≥' + str(labels[-1])
    plt.gca().set_xticklabels(labels)

    box_pairs = [(labels[c], labels[c+1]) for c, _ in enumerate(labels[:-1])]
    add_stat_annotation_boxplot_no_hue(plt.gca(), data, xcol='quantiles', ycol='name', 
                                        order=labels,
                                        box_pairs=box_pairs,
                                        ymax=12, delta=.05, h=.1,
                                    )
    plt.title("H3K27ac in metaloop anchors")
    fig.savefig('./plots/paper/s28/h3k27ac_metaloops.pdf', bbox_inches = 'tight')

def fraction_tss_with_metaloops(my_tss_df, metaloop_anchors, gene_dict):
    l = pbt.BedTool.from_dataframe(my_tss_df.iloc[:, :4]).slop(b=2_500, genome='mm10').intersect(metaloop_anchors, c=True
                                                                        ).to_dataframe()
    rpkm_bins = pd.qcut(gene_dict['Resting']['rpkm'].sort_values(),
            np.linspace(0, 1, 10), labels=False)
    rpkm_bins = pd.qcut(gene_dict['Resting']['rpkm'].sort_values(),
            np.linspace(0, 1, 10), labels=False).dropna()
    xs = np.unique(rpkm_bins)
    frac_ys = []
    ys = []
    for u in xs:
        genes = rpkm_bins.index[rpkm_bins == u]
        values = l.set_index('name').loc[genes]['score']
        y = np.nanmean(values > 0)
        ys.append(y)
    ys = arr(ys)
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)
    scatter_with_pearson(xs, ys, axs, marker='o')
    plt.plot(xs, ys)

    plt.xlabel("TSS quantiles")
    plt.ylabel("Fraction with metaloops")
    add_xaxis_labels("Low Expr", 'High Expr', plt.gca(), fontsize=6)
    plt.gca().set_xticklabels([])
    plt.ylim([0, .6])
    plt.title("Fraction TSS with metaloops")
    fig.savefig('./plots/paper/s28/tss_with_metaloops.pdf', bbox_inches = 'tight')    
    
from plot_pvals import add_stat_annotation_boxplot_no_hue
def rpkm_by_metaloops(metaloop_anchors, my_tss_df, gene_dict):
    l = pbt.BedTool.from_dataframe(my_tss_df.iloc[:, :4]
                                ).slop(b=2_500, genome='mm10').saveas().intersect(metaloop_anchors, c=True).to_dataframe(header=None)
    l['quantiles'] = pd.cut(l['score'], [-np.inf, .5, 1.5, 2.5,np.inf], labels=False)
    xs = []
    ys = [] 
    ns = []
    colors = sns.color_palette("coolwarm", n_colors = 5)
    # for c, (clow, chigh) in enumerate(cutoffs):
    # names = l[3][(l['score'] >= clow) & (l['score'] <= chigh)]
    x = l.set_index('name')['quantiles']
    y = (gene_dict['Resting']['rpkm'].loc[x.index])
    data = pd.concat([x, y], axis=1)
    fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)

    sns.boxplot(data=data, x='quantiles', y='rpkm', showfliers=False, ax=plt.gca())
    add_stat_annotation_boxplot_no_hue(plt.gca(), data, xcol='quantiles', ycol='rpkm', 
                                    order=[0, 1, 2, '≥3'], 
                                    box_pairs=[[0, 1], [1, 2], [2, '≥3'], ],
                                    ymax=1.6, delta=.05, h=.02,
                                    )

    plt.gca().set_xticklabels([0, 1, 2, '≥3'])
    plt.ylabel("RPKM")
    plt.xlabel("# Metaloops")
    fig.savefig('./plots/paper/s28/rpkm_in_metaloops.pdf', bbox_inches = 'tight')
    plt.title("RPKM by metaloops")    
# def atac_peaks_in_metaloop_anchors(metaloop_anchors, atac_peaks):
#     metaloops_with_count = pbt.BedTool.from_dataframe(metaloop_anchors.to_dataframe().value_counts().reset_index())
#     l2 = add_chr_to_bedtool(metaloops_with_count
#                             ).intersect(add_chr_to_bedtool(atac_peaks), c=True).to_dataframe(header=None)

#     l2['quantiles'] = l2['name']
#     x = l2['quantiles']
#     y = l2['score']
#     data = pd.concat([x, y], axis=1)

#     data['score'] = data['score'].clip(0, 2)
#     data['quantiles'] = data['quantiles'].clip(0, 4)
#     s = data.value_counts(['quantiles', 'score']).unstack().fillna(0)
#     s = s.div(s.sum(axis=1), axis=0)

#     fig, axs = init_subplots_exact(1, 1, fgsz=(40*mm, 40*mm), dpi = 100)

#     s.plot.bar(stacked=True, ax=plt.gca())
#     plt.xlabel("# Metaloops")
#     plt.ylabel("# peaks")

#     box_pairs = [(labels[c], labels[c+1]) for c, _ in enumerate(labels[:-1])]

#     labels = [int(x.get_text()) for x in plt.gca().get_xticklabels()]
#     labels[-1] = '≥' + str(labels[-1])
#     plt.gca().set_xticklabels(labels)

#     box_pairs = [(labels[c], labels[c+1]) for c, _ in enumerate(labels[:-1])]
#     add_stat_annotation_boxplot_no_hue(plt.gca(), data, xcol='quantiles', ycol='score', 
#                                         order=labels,
#                                         box_pairs=box_pairs,
#                                         ymax=1, delta=.05, h=.02,
#                                     )

#     plt.xticks(rotation=0)
#     plt.legend(bbox_to_anchor=(1, 1))
#     plt.title("# ATAC peaks in metaloop anchors")
#     fig.savefig('./plots/paper/s28/atac_peaks_in_metaloops.pdf', bbox_inches = 'tight')    