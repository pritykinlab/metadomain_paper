import statsmodels
import statsmodels.stats
import statsmodels.stats.multitest
from plotting_functions import * 

good_foxp3_experiments = [
 '_aTreg_Foxp3_ChIP_Rep2__Mus_musculus__ChIP-Seq',
 '_aTreg_Foxp3_ChIP_Rep3__Mus_musculus__ChIP-Seq',
 '_Foxp3_ChIP-seq_bulk_Treg_rep1__Mus_musculus__ChIP-Seq',
 '_Foxp3_ChIP-seq_in_Foxp1+_Treg__rep2__Mus_musculus__ChIP-Seq',
 '_Foxp3_ChIP-seq_bulk_Treg_rep2__Mus_musculus__ChIP-Seq',
 '_Foxp3_ChIP-seq_in_Foxp1+_Treg__rep1__Mus_musculus__ChIP-Seq',
 '_aTreg_Foxp3_ChIP_Rep4__Mus_musculus__ChIP-Seq',
 '_aTreg_Foxp3_ChIP_Rep1__Tech_Rep_1___Mus_musculus__ChIP-Seq',
 '_aTreg_Foxp3_ChIP_Rep1__Tech_Rep_2___Mus_musculus__ChIP-Seq',
 '_Treg_Foxp3_ChIP_Rep2__Mus_musculus__ChIP-Seq',
 '_Treg_Foxp3_ChIP_Rep1__Tech_Rep1___Mus_musculus__ChIP-Seq',
 '_Treg_Foxp3_ChIP_Rep1__Tech_Rep2___Mus_musculus__ChIP-Seq',
 'Foxp3_ChIP-seq_of_Treg',
'_Foxp3_ChIP-seq_bulk_Treg_rep3__Mus_musculus__ChIP-Seq'
]
import pandas as pd
import scipy

def get_cellstatus_from_experiments(index):
    celltypes = pd.Series(index).copy()
    celltypes.index = index
    statusdict = {
        'a' : (index.str.contains('activated')) | (index.str.contains('aT')) ,
        'r' : (index.str.contains('rest')) | (index.str.contains('rT')) ,
        't' : (index.str.contains('tT')),
        'm' : (index.str.contains('mT')),
    }
    for i, v in statusdict.items():
        celltypes[v] = i
    celltypes[celltypes==index] = 'bulk'
    return celltypes

def get_celltype_from_experiments(index):
    celltypes = pd.Series(index).copy()
    celltypes[celltypes.index] = 'Not_Treg_Tcon'
    celltypes.index = index
    statusdict = {
        'CD4' : ((index.str.contains('CD4_')) & (~(index.str.contains('CD4SP')))) | (index.str.contains('_H3K9me3_WT__Mus_musculus__ChIP')),
        'CD4SP' : (index.str.contains('CD4SP')),

        'Treg Prec.' : (index.str.contains('precursor')),
        'Treg' : (index.str.contains('Tr')) & (~(index.str.contains('precursor'))),
        'Tcon' : (index.str.contains('Tcon')) | (index.str.contains('conv')),

    }
    for i, v in statusdict.items():
        celltypes[v] = i
    outliers = [
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep2__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep1__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep3__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep2__Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_Foxp1+_Treg__rep1__Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_Foxp1-_Treg__genetic_control___Mus_musculus__ChIP-Seq',
               ]
    celltypes[celltypes.index.isin(outliers)] = "outlier"
    celltypes[celltypes.index.str.contains("KO")] = 'KO'
    celltypes[celltypes.index.str.contains("CNS")] = 'CNSKO'
    celltypes[celltypes == '_Tn_H3K27me3_ChIP_Rep1__Mus_musculus__ChIP-Seq'] = 'Tcon'
    celltypes[celltypes == '_Tn_H3K27me3_ChIP_Rep2__Mus_musculus__ChIP-Seq'] = 'Tcon'
    celltypes[celltypes.index.str.contains("DN2")] = 'DN2'
    print(1)
    celltypes.name = 'Celltype'
    return celltypes

def get_target_from_experiments(index):
    inddict = {
        'Stat5' : index.str.contains("Stat5"),
        'Ets1' : index.str.contains("Ets1"),
        'Tet2' : index.str.contains("Tet2"),
        'CREB' : index.str.contains("CREB"),
        'Elf' : index.str.contains("Elf"),
        'H3K4me1' : index.str.contains("H3K4me1"),
        'H3K4me3' : index.str.contains("H3K4me3"),
        'H3K27ac' : index.str.contains("H3K27ac") | index.str.contains("H3K27Ac"),
        'H3K27me3' : index.str.contains("H3K27me3"),
        'Smc1' : index.str.contains("Smc1"),
        'Foxp3' : index.isin(good_foxp3_experiments) | index.str.contains("Foxp3_CUTRUN"),
        'Foxp1' : index.str.contains("Foxp1_ChI"),
        'Ctcf' : index.str.contains("CTCF"),        
        'TSS' : index.str.contains("tss"),        
        'Relap65' : index.str.contains("Relap65"),
        'Runx1' : index.str.contains("Runx1"),        
        'H3K9me3' : index.str.contains("H3K9me3"),        
        'Cbfb' : index.str.contains("Cbfb"),
        'Ep300' : index.str.contains("Ep300"),
        'Input' : index.str.contains("Control") | index.str.contains("Input") | index.str.contains("control"),
        'Bcl11b' : index.str.contains("Bcl11b"),
        'Med1' : index.str.contains("Med1"),
        'Satb1' : index.str.contains("Satb1"),
        'TCF1' : index.str.contains("TCF1-CUT"),
        'Lef1' : index.str.contains("Lef1"),
        'cJun' : index.str.contains("cJun"),
        'IRF4' : index.str.contains("IRF4"),
        'ATAC' : index.str.contains("ATAC"),
    }        
    experiments = pd.DataFrame([], index=index).copy()
    experiments['Experiment'] = "None"
    for i, protein in inddict.items():
        experiments['Experiment'].iloc[protein] = i
    outliers = [
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep2__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep1__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep3__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep2__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep1__Tech_Rep_1___Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep1__Tech_Rep_2___Mus_musculus__ChIP-Seq'
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep2__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Treg_Foxp3_ChIP_Rep1__Tech_Rep1___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep2__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep2__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep1__Tech_Rep_1___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep2__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Treg_Foxp3_ChIP_Rep1__Tech_Rep1___Mus_musculus__ChIP-Seq',
                '_Treg_Foxp3_ChIP_Rep1__Tech_Rep2___Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep1__Tech_Rep_2___Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_conventional_T_cells__rep2__Mus_musculus__ChIP-Seq-1',
                '_CD4_Cbfb_ChIP__Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_Foxp1-_Treg__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_Foxp1+_Treg__rep2__Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep3__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_Foxp1+_Treg__rep1__Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_Foxp1+_Treg__rep3__Mus_musculus__ChIP-Seq',
                '_Foxp1_ChIP-seq_in_conventional_T_cells__rep2__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep2__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep2__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_aTreg_Foxp3_ChIP_Rep1__Tech_Rep_1___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep2__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Treg_Foxp3_ChIP_Rep1__Tech_Rep1___Mus_musculus__ChIP-Seq',
               ]
    experiments['Experiment'][experiments.index.isin(outliers)] = 'No'
    experiments['Experiment'].replace('None', 'No', inplace=True)
    return experiments


def adata_calculate_differences_in_enhancer_subsets(adata, indsoi, indsoi2,
                                                             uplabel='Higher in \nmegaloop anchors',
                                                             downlabel='Higher in \nother anchors',
                                                             zscore = False):
    inddict = {
        'Stat5' : adata.var.index.str.contains("Stat5"),
        'Ets1' : adata.var.index.str.contains("Ets1"),
        'Tet2' : adata.var.index.str.contains("Tet2"),
        'CREB' : adata.var.index.str.contains("CREB"),
        'Elf' : adata.var.index.str.contains("Elf"),
        'H3K4me1' : adata.var.index.str.contains("H3K4me1"),
        'H3K4me3' : adata.var.index.str.contains("H3K4me3"),
        'H3K27ac' : adata.var.index.str.contains("H3K27ac") | adata.var.index.str.contains("H3K27Ac"),
        'H3K27me3' : adata.var.index.str.contains("H3K27me3"),
        'Smc1' : adata.var.index.str.contains("Smc1"),
        'Foxp3' : adata.var.index.isin(good_foxp3_experiments) | adata.var.index.str.contains("Foxp3_CUTRUN"),
        'Foxp1' : adata.var.index.str.contains("Foxp1_ChI"),
        'Ctcf' : adata.var.index.str.contains("CTCF"),        
        'TSS' : adata.var.index.str.contains("tss"),        
        'Relap65' : adata.var.index.str.contains("Relap65"),
        'Runx1' : adata.var.index.str.contains("Runx1"),        
        'H3K9me3' : adata.var.index.str.contains("H3K9me3"),        
        'Cbfb' : adata.var.index.str.contains("Cbfb"),
        'Ep300' : adata.var.index.str.contains("Ep300"),
        'Input' : adata.var.index.str.contains("Control") | adata.var.index.str.contains("Input") | adata.var.index.str.contains("control"),
        'Bcl11b' : adata.var.index.str.contains("Bcl11b"),
        'Med1' : adata.var.index.str.contains("Med1"),
        'Satb1' : adata.var.index.str.contains("Satb1"),
        'TCF1' : adata.var.index.str.contains("TCF1-CUT"),
        'Lef1' : adata.var.index.str.contains("Lef1"),
        'cJun' : adata.var.index.str.contains("cJun"),
        'IRF4' : adata.var.index.str.contains("IRF4"),
        'ATAC' : adata.var.index.str.contains("ATAC"),
    }    
    subadata = adata[indsoi | indsoi2, :]
    if zscore:
        subadata.X = scipy.stats.zscore(subadata.X, axis=0)

    delta_df = pd.DataFrame()
    delta_df['chip_mega'] = subadata[adata.obs.index[indsoi], :].X.mean(axis=0)
    delta_df['chip_no_mega'] = subadata[adata.obs.index[indsoi2], :].X.mean(axis=0)
    delta_df['delta'] = delta_df['chip_mega'] - delta_df['chip_no_mega']
    delta_df.index = subadata.var.index

    metadata_df = pd.DataFrame([], index=subadata.var.index)
    val_df = pd.DataFrame
    delta_df = pd.DataFrame()
    m1 = subadata[adata.obs.index[indsoi], :].X.copy()
    m2 = subadata[adata.obs.index[indsoi2], :].X.copy()
    # baseline = adata.X.mean(axis=0)
    # basevar = adata.X.std(axis=0)
    delta_df['chip_mega'] = m1.mean(axis=0)
    delta_df['chip_no_mega'] = m2.mean(axis=0)
    pvalues = scipy.stats.mannwhitneyu(m1, m2)[1]
    delta_df['delta'] = (delta_df['chip_mega'] - delta_df['chip_no_mega'])#/basevar
    delta_df['pvalue'] = pvalues
    delta_df.index = subadata.var.index
    sorted_vals = (delta_df['delta']).sort_values()[::-1]
    delta_df['Experiment'] = "None"
    for i, protein in inddict.items():
        delta_df.loc[protein, 'Experiment'] = i
    outliers = ['H3K27ac_ChIP-seq_of_ThpokCreSatb1CKO_Foxp3+Tconv',
               'H3K27ac_ChIP-seq_of_ThpokCreSatb1CKO_Foxp3-Tconv',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep2__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep1__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep3__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep2__Mus_musculus__ChIP-Seq',
                '901R3b_Stat5ChIP_Treg_periphery_CNS0KO_rep3',
                '901R3d_Stat5ChIP_Treg_periphery_CNS03dKO_rep3'
               ]
    outliers += list(delta_df.index[delta_df.index.str.contains("CNS")])
    delta_df['Experiment'][delta_df.index.isin(outliers)] = 'No'
    delta_df['Experiment'].replace('None', 'No', inplace=True)
    delta_df = delta_df[~(delta_df.Experiment=='No')]
    delta_df['Celltype'] = 'None'
    return delta_df, (adata[indsoi, :], adata[indsoi2, :])


import statsmodels
import statsmodels.stats
import statsmodels.stats.multitest
import seaborn as sns
import matplotlib.pyplot as plt
def adata_calculate_differences_in_enhancer_subsets_stripplot(adata, indsoi, indsoi2,
                                                             uplabel='Higher in \nmegaloop anchors',
                                                             downlabel='Higher in \nother anchors',
                                                             zscore=False
                                                             ):
    subadata = adata[indsoi | indsoi2, :]
    if zscore:
        subadata.X = scipy.stats.zscore(subadata.X, axis=0)

    delta_df = pd.DataFrame()
    delta_df['chip_mega'] = subadata[indsoi, :].X.mean(axis=0)
    delta_df['chip_no_mega'] = subadata[indsoi2, :].X.mean(axis=0)
    delta_df['delta'] = delta_df['chip_mega'] - delta_df['chip_no_mega']
    delta_df.index = subadata.var.index
    sorted_vals = (delta_df['delta']).sort_values()[::-1]
    inddict = {
        'Stat5' : sorted_vals.index.str.contains("Stat5"),
        'Ets1' : sorted_vals.index.str.contains("Ets1"),
        'Tet2' : sorted_vals.index.str.contains("Tet2"),
        'CREB' : sorted_vals.index.str.contains("CREB"),
        'Elf' : sorted_vals.index.str.contains("Elf"),
        'H3K4me1' : sorted_vals.index.str.contains("H3K4me1"),
        'H3K4me3' : sorted_vals.index.str.contains("H3K4me3"),
        'H3K27ac' : sorted_vals.index.str.contains("H3K27ac") | sorted_vals.index.str.contains("H3K27Ac"),
        'H3K27me3' : sorted_vals.index.str.contains("H3K27me3"),
        'Smc1' : sorted_vals.index.str.contains("Smc1"),
        'Foxp3' : sorted_vals.index.isin(good_foxp3_experiments) | sorted_vals.index.str.contains("Foxp3_CUTRUN"),
        'Foxp1' : sorted_vals.index.str.contains("Foxp1_ChI"),
        'Ctcf' : sorted_vals.index.str.contains("CTCF"),        
        'TSS' : sorted_vals.index.str.contains("tss"),        
        'Relap65' : sorted_vals.index.str.contains("Relap65"),
        'Runx1' : sorted_vals.index.str.contains("Runx1"),        
        'H3K9me3' : sorted_vals.index.str.contains("H3K9me3"),        
        'Cbfb' : sorted_vals.index.str.contains("Cbfb"),
        'Ep300' : sorted_vals.index.str.contains("Ep300"),
        'Input' : sorted_vals.index.str.contains("Control") | sorted_vals.index.str.contains("Input") | sorted_vals.index.str.contains("control"),
        'Bcl11b' : sorted_vals.index.str.contains("Bcl11b"),
        'Med1' : sorted_vals.index.str.contains("Med1"),
        'Satb1' : sorted_vals.index.str.contains("Satb1"),
        'TCF1' : sorted_vals.index.str.contains("TCF1-CUT"),
        'Lef1' : sorted_vals.index.str.contains("Lef1"),
        'cJun' : sorted_vals.index.str.contains("cJun"),
        'IRF4' : sorted_vals.index.str.contains("IRF4"),
    }
    df = pd.DataFrame(sorted_vals)
    df['Experiment'] = "None"
    for i, protein in inddict.items():
        df.loc[protein, 'Experiment'] = i
    outliers = ['H3K27ac_ChIP-seq_of_ThpokCreSatb1CKO_Foxp3+Tconv',
               'H3K27ac_ChIP-seq_of_ThpokCreSatb1CKO_Foxp3-Tconv',
                '_Foxp3_ChIP-seq_in_Foxp1+_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep1__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep3__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_conventional_T_cells__rep2__genetic_control___Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep1__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep3__Mus_musculus__ChIP-Seq',
                '_Foxp3_ChIP-seq_in_Foxp1-_Treg__rep2__Mus_musculus__ChIP-Seq',
               '901R3b_Stat5ChIP_Treg_periphery_CNS0KO_rep3',
                '901R3d_Stat5ChIP_Treg_periphery_CNS03dKO_rep3'
               ]
    outliers += list(df.index[df.index.str.contains("CNS")])
    df.loc[df.index.isin(outliers), 'Experiment'] = 'No'
    df['Experiment'].replace('None', 'No', inplace=True)
    df = df[~(df.Experiment=='No')]
    df.index=df['Experiment'].values
    df = df.loc[df.groupby("Experiment").median().sort_values("delta", ascending=False).index]

    fig, axs = init_subplots_exact(1, 1, fgsz=(80*mm, 40*mm), dpi=300)
    ax = axs
    sns.boxplot(data=df, x = 'Experiment', y = 'delta', ax=ax, showfliers=False)
    sns.stripplot(data=df, x = 'Experiment', y = 'delta', ax=ax, color='black', s=2)
    plt.xticks(rotation=90);
    plt.ylim([-1.2, 1.2])


    ax.set_title("Differential ChIP signal")
    ax.set_xlabel("Factor")
    ax.text(-.1, 1.1, uplabel, transform=ax.transAxes, ha='center', va='center', fontsize=6)
    ax.text(-.1, -.1, downlabel, transform=ax.transAxes, ha='center', va='center', fontsize=6)
    return fig, df, delta_df