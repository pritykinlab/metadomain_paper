import glob
import pybedtools as pbt
import sys
from aux_functions import get_col

def init_gene_dicts(path='./peaks/*thresh*'):
    gene_dict = {}
    for file in glob.glob(path):
        cond1, cond2 = file.split("/")[-1].split("_thresh")[0].split("_vs_")
        thresh = file.split("_thresh=")[1].split(".csv")[0]
        bedtool = pbt.BedTool(file)
        gene_dict[f'{cond1}_vs_{cond2}_thresh={thresh}'] = bedtool

    gene_set_dict = {}
    for cond, bedfile in gene_dict.items():
        remainder, thresh = cond.split("_thresh=")
        cond1, cond2 = remainder.split("_vs_")
        cond1_up = bedfile.filter(lambda x: (float(x[4]) < .05) & (float(x[6]) > 0)).saveas()
        cond2_up = bedfile.filter(lambda x: (float(x[4]) < .05) & (float(x[6]) < 0)).saveas()
        gene_set_dict[f'{cond1}_{cond}'] = get_col(cond1_up, 3)
        gene_set_dict[f'{cond2}_{cond}'] = get_col(cond2_up, 3)
    return gene_dict, gene_set_dict

def init_lfc_dicts(adata):
    gene_lfc_dict = {}
    gene_name_dict = {}
    gene_pval_dict = {}
    for cond, bedfile in gene_dict.items():
        names = get_col(bedfile, 3)
        lfcs = get_col(bedfile, 6).astype(float)
        pvals = get_col(bedfile, 4).astype(float)
        gene_lfc_dict[cond] = lfcs[[x in adata.var.index for x in names]]
        gene_name_dict[cond] = names[[x in adata.var.index for x in names]]
        gene_pval_dict[cond] = pvals[[x in adata.var.index for x in names]]
    return gene_lfc_dict, gene_name_dict, gene_pval_dict

import pandas as pd
def load_lfc_and_pval_df(globpath):
    lfc_df = []
    wald_df = []
    pval_df = []
    names = []
    basemeans = []
    for file in glob.glob(globpath):
        name = file.split('/')[-1].split("_thresh")[0]
        df = pd.read_csv(file, sep=' ')
        lfc_df.append(df['log2FoldChange'])
        pval_df.append(df['padj'])
        basemeans.append(df['baseMean'])
        wald_df.append(df['log2FoldChange'] / df['lfcSE'])
        names.append(name)
    lfc_df = pd.DataFrame(lfc_df).T
    pval_df = pd.DataFrame(pval_df).T
    wald_df = pd.DataFrame(wald_df).T

    lfc_df.columns = names
    wald_df.columns = names
    pval_df.columns = names
    basemean_df = pd.DataFrame(basemeans).T
    basemean_df.columns = names
    return lfc_df, pval_df, basemean_df, wald_df