import pandas as pd
from aux_functions import make_int, index_to_bedtool, bedtool_to_index
import pybedtools as pbt
def load_deseq_df(path, drop_overlaps=False, **kwargs):
    df = pd.read_csv(path, sep=' ', **kwargs)
    df = mod_index(df)
    if drop_overlaps:
        rows = list(index_to_bedtool(df.index).sort())
        n_overlaps = 1
        n_loops = 0
        while n_overlaps > 0:
            n_overlaps = 0
            for c in range(len(rows)-1):
                if c >= len(rows)-1:
                    continue
                row = make_int(rows[c])
                nextrow = make_int(rows[c+1])
                if (row[0] == nextrow[0]) and (row[2] > nextrow[1]):     
                    del rows[c+1]
                    n_overlaps += 1
            print(f"Dropped {n_overlaps}")
            n_loops += 1
        return df.loc[rows]

    return df

def drop_overlaps(index):
    rows = list(index_to_bedtool(index).sort())
    n_overlaps = 1
    n_loops = 0
    while n_overlaps > 0:
        n_overlaps = 0
        for c in range(len(rows)-1):
            if c >= len(rows)-1:
                continue
            row = make_int(rows[c])
            nextrow = make_int(rows[c+1])
            if (row[0] == nextrow[0]) and (row[2] > nextrow[1]):     
                del rows[c+1]
                n_overlaps += 1
        print(f"Dropped {n_overlaps}")
        n_loops += 1
    return bedtool_to_index(rows)

def mod_index(df):
    df = df.copy()
    df.index = ["_".join([make_int(x.split("_"))[0], str(make_int(x.split("_"))[1]+1), str(make_int(x.split("_"))[2])])
            for x in df.index ]
    return df

def unmod_index(df):
    df = df.copy()
    df.index = ["_".join([make_int(x.split("_"))[0], str(make_int(x.split("_"))[1]-1), str(make_int(x.split("_"))[2])])
            for x in df.index ]
    return df


import glob
def parse_DESeq2_dfs(dfdict):
    rows = []
    lfc_df = pd.DataFrame()
    pval_df = pd.DataFrame()
    basemean_df = pd.DataFrame()
    wald_df = pd.DataFrame()

    lfc_cols = []
    pval_cols = []
    basemean_cols = []
    wald_df_cols = []
    cols = []
    for cond, df in dfdict.items():
        lfc_cols.append(df['log2FoldChange'])
        pval_cols.append(df['padj'])
        basemean_cols.append(df['baseMean'])
        wald_df_cols.append(df['log2FoldChange'] / df['lfcSE'])
        cols.append(cond)
        
    lfc_df = pd.concat(lfc_cols, axis=1)
    pval_df = pd.concat(pval_cols, axis=1)
    basemean_df = pd.concat(basemean_cols, axis=1)
    wald_df = pd.concat(wald_df_cols, axis=1)
    
    lfc_df.columns = cols
    pval_df.columns = cols
    basemean_df.columns = cols
    wald_df.columns = cols    
    
    differential_gene_counts = pd.DataFrame(rows, columns = ['cond1', 'cond2', 'name', 'cond', 'down', 'ns', 'up'])    
    return lfc_df, pval_df, basemean_df, wald_df, differential_gene_counts

