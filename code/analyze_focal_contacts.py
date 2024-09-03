import scipy
import scipy.signal
import numpy as np
import pandas as pd
from hic_zscore_functions import make_obs_exp_nolog
from copy import deepcopy
arr = np.asarray


def make_treg_tcon_df(indsoi, valdf, treg_specific, tcon_specific, int_to_chip_clust):
    df_treg = pd.DataFrame(); df_tcon = pd.DataFrame()
    keys = []
    for i in indsoi:        
        treg_specific_loops = np.where(treg_specific[i])[0]
        tcon_specific_loops = np.where(tcon_specific[i])[0]
        treg_subdf = valdf[valdf['Inds'].isin(treg_specific_loops)]
        tcon_subdf = valdf[valdf['Inds'].isin(tcon_specific_loops)]

        treg_clusts = treg_subdf['ChIP_Clust']
        tcon_clusts = tcon_subdf['ChIP_Clust']; n = 3;

        vs, counts = np.unique(treg_clusts, return_counts=True);
        treg_counts = np.zeros(n)
        treg_counts[vs] = counts
        treg_counts = treg_counts

        vs, counts = np.unique(tcon_clusts, return_counts=True)
        tcon_counts = np.zeros(n)
        tcon_counts[vs] = counts
        tcon_counts = tcon_counts

        keys = [int_to_chip_clust[x] for x in range(n)]

        d_treg = (treg_counts)#/treg_counts.sum()
        d_tcon = (tcon_counts)#/tcon_counts.sum()

        ind_subdf = valdf[valdf['Inds'].isin([i])]
        h3k27ac_delta = np.log(ind_subdf['h3k27ac_treg']/ind_subdf['h3k27ac_tcon']).values
        h3k27me3_delta = np.log(ind_subdf['h3k27me3_treg']/ind_subdf['h3k27me3_tcon']).values

        cluster_delta = d_treg
        row = pd.DataFrame([list(h3k27ac_delta) + list(h3k27me3_delta) + list(d_treg)])
        row.index = [i]
        df_treg = pd.concat([df_treg, row], axis=0)

        cluster_delta = d_treg
        row = pd.DataFrame([list(h3k27ac_delta) + list(h3k27me3_delta) + list(d_tcon)])
        row.index = [i]
        df_tcon = pd.concat([df_tcon, row], axis=0)

    df_treg.columns = ['∆ H3K27ac', '∆ H3K27me3' ] + list(keys)
    df_tcon.columns = ['∆ H3K27ac', '∆ H3K27me3' ] + list(keys)
    return df_treg, df_tcon




