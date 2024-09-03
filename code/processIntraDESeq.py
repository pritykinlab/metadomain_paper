import numpy as np
import pandas as pd

def process_deseq_df_to_mats(df, all_ind_to_region, pcolumn = 'padj'):
    n = len(all_ind_to_region)
    deseq_effect_mat = np.zeros((n, n))*np.nan
    deseq_pval_mat = np.zeros((n, n))*np.nan
    deseq_lfc_mat = np.zeros((n, n))*np.nan
    lfc = df['log2FoldChange']
    effect = df['log2FoldChange'] / df['lfcSE']
    pvals = df[pcolumn]

    arr = np.asarray
    tmp = df.reset_index()['index'].apply(lambda x: x.split('_')).values
    rows, cols = arr(list(zip(*tmp)))
    rows, cols = rows.astype(int), cols.astype(int)
    deseq_lfc_mat[rows, cols] = lfc
    deseq_lfc_mat[cols, rows] = lfc

    
    deseq_effect_mat[rows, cols] = effect
    deseq_effect_mat[cols, rows] = effect

    deseq_pval_mat[rows, cols] = pvals
    deseq_pval_mat[cols, rows] = pvals    
    return deseq_lfc_mat, deseq_effect_mat, deseq_pval_mat
