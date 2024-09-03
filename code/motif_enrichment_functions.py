import scipy

import scipy.ndimage
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pybedtools as pbt
from aux_functions import fetch_mean_matched_values

def get_motif_enrichments(df, motif_df, mean_match=True, mod=10):
    rs = []
    for c, col in enumerate(motif_df):
        motif_idx = motif_df[col] > 0

        x = df.loc[motif_idx, 'log2FoldChange']
        if mean_match:
            motif_nonidx = fetch_mean_matched_values(df.loc[motif_idx, 'baseMean'],
                                       df.loc[~motif_idx, 'baseMean']).index
        else:
            motif_nonidx = motif_df[col] == 0
        y = df.loc[motif_nonidx, 'log2FoldChange']
        _, p = scipy.stats.ranksums(x, y)
        delta = (x.mean() - y.mean())
        rs.append([col, p, delta])    
        if c % mod == 0:
            print("Done with", c)
    return pd.DataFrame(rs, columns = ['factor', 'p', 'delta'])


from concurrent.futures import ProcessPoolExecutor

# Ensure your fetch_mean_matched_values function is defined here

def compute_enrichment(args):
    """
    Compute enrichment for a single motif.
    This function is designed to be run in parallel.
    """
    df, motif_df_path, col = args
    motif_df = pd.read_csv(motif_df_path, index_col = 0)
    motif_df = motif_df.loc[df.index]
    motif_df = motif_df[~motif_df.index.duplicated()]
    motif_idx = motif_df[col] > 0
    x = df.loc[motif_idx, 'log2FoldChange']
    if motif_idx.mean() > .7:
        return [col, np.nan, np.nan]
    elif motif_idx.sum() > .3:
        motif_nonidx = motif_df[col] == 0
    else:
        motif_nonidx = fetch_mean_matched_values(df.loc[motif_idx, 'baseMean'], df.loc[~motif_idx, 'baseMean']).index
    y = df.loc[motif_nonidx, 'log2FoldChange']
    _, p = scipy.stats.ranksums(x, y)
    delta = (x.mean() - y.mean())
    print("Done with", col)
    return col, p, delta

def get_motif_enrichments_parallel(df, motif_df_path, n_cores=20):
    motif_df = pd.read_csv(motif_df_path, index_col=0)
    # Prepare arguments for parallel computation
    tasks = [(df, motif_df_path, col) for col in motif_df.columns]

    # Process in parallel
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = list(executor.map(compute_enrichment, tasks))

    # Return results as a DataFrame
    return pd.DataFrame(results, columns=['factor', 'p', 'delta'])

def get_motif_enrichments_using_idx(motif_df, idx1, idx2):
    rows = []
    for col in motif_df:
        x = motif_df.loc[idx1, col]
        y = motif_df.loc[idx2, col]
        _, p = scipy.stats.fisher_exact([[(x>0).sum(), (x==0).sum()], [(y>0).sum(), (y==0).sum()]])
        xval = (x.sum())
        yval = (y.sum())
        if (xval==0):
            xval = 1
        if (yval==0):
            yval = 1
        xval = xval / len(x)
        yval = yval / len(y)
        delta = np.log2(xval/yval)
        rows.append([col, p, delta])
    return pd.DataFrame(rows, columns=['factor', 'p', 'delta'])






