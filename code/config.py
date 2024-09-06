row_colors = ['lightgreen', 'green', 'orange']

columns_to_names = {
    0 : 'Constitutive',
    4 : 'Dynamic',
    18 : 'Repressive',
}

row_colors_dict = {
'Constitutive' : 'lightgreen',
'Dynamic' : 'green',
'Repressive' : 'orange',
}


from aux_functions import bedtool_to_index, add_chr_to_bedtool
import pybedtools as pbt
import pandas as pd
def get_h3k27ac_motif_counts(pco=5e-5):    
    atac_peak_motif_counts = bedtool_to_index(add_chr_to_bedtool(pbt.BedTool('h3k27ac/all_threshold_27ac.csv',
                                )))
    atac_peak_motif_counts.index = atac_peak_motif_counts.values
    meme_motif_df = pd.read_csv('h3k27ac_peaks/h3k27ac_motifs.bed', sep='\t')
    motif_id_to_name_dict = dict(zip(meme_motif_df['motif_id'], meme_motif_df['dedup_protein_name']))
    meme_motif_df = meme_motif_df[meme_motif_df['pval'] < pco]
    meme_motif_df = meme_motif_df.value_counts(['peak', 'dedup_protein_name']).unstack()
    meme_motif_df.index = meme_motif_df.index.str.lower().str.replace(":", "_").str.replace("-", "_")
    meme_motif_df = meme_motif_df.fillna(0)
    vs = []
    missing_inds = atac_peak_motif_counts.index[~atac_peak_motif_counts.index.isin(meme_motif_df.index)]
    if len(missing_inds) > 0:
        for _ in missing_inds:
            v = pd.Series(index=meme_motif_df.columns).fillna(0)
            vs.append(v)
        
        tmpdf = pd.concat(vs, axis=1).T
        tmpdf.index = missing_inds
        meme_motif_df = pd.concat([meme_motif_df, tmpdf], axis=0)
    meme_motif_df = meme_motif_df.loc[atac_peak_motif_counts.index]
    return meme_motif_df, motif_id_to_name_dict


import numpy as np
def load_scrna():
    susie_path = './coexpression_analysis/'
    # treg_precursor_rna_umi = pd.read_parquet(susie_path + 'GSM3978654_CD4plusCD25plus_dense.umis_per_gene.parquet')
    # tregs_rna_umi = pd.read_parquet(susie_path + 'GSM3978655_Foxp3plus_dense.umis_per_gene.parquet')

    treg_precursor_rna_corr = pd.read_parquet('./coexpression_analysis/GSM3978654_CD4plusCD25plus_dense.coexpression_mat.parquet')
    tregs_rna_corr = pd.read_parquet('./coexpression_analysis/GSM3978655_Foxp3plus_dense.coexpression_mat.parquet')

    treg_precursor_rna_corr += np.diag(np.diag(treg_precursor_rna_corr)*np.nan)
    tregs_rna_corr += np.diag(np.diag(tregs_rna_corr)*np.nan)
    same_genes = treg_precursor_rna_corr.index.intersection(tregs_rna_corr.index)

    df1 = pd.read_parquet('./coexpression_analysis/GSM3978654_CD4plusCD25plus_dense.umis_per_gene.parquet')
    df2 = pd.read_parquet('./coexpression_analysis/GSM3978655_Foxp3plus_dense.umis_per_gene.parquet')

    df1 = df1.loc[same_genes]
    df2 = df2.loc[same_genes]

    same_genes = same_genes[(df1['Raw_UMI_count'] > 50) & (df2['Raw_UMI_count'] > 50)]


    treg_precursor_rna_corr = treg_precursor_rna_corr.loc[same_genes, same_genes]
    tregs_rna_corr = tregs_rna_corr.loc[same_genes, same_genes]
    return treg_precursor_rna_corr, tregs_rna_corr 

cp /Genomics/argo/users/gdolsten/pritlab/yanjin_motifreg_hic/meme/fimo_output/treg_tcon_atac.bed atac_peaks/atac_motifs.bed


from load_DESeq2_data import mod_index
def get_motif_counts(atac_peaks_index, pco=5e-5):
    # atac_peak_motif_counts = pd.read_csv('/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/rudensky_scrna/prelim-analysis/call_motifs/output/counts/previous_bulk_atac_peaks_Mus_musculus/previous_bulk_atac_peaks_p=1e-05.csv',
                                #   index_col = 0) > 0 
    meme_motif_df = pd.read_csv('./atac/atac_motifs.bed', sep='\t')
    motif_id_to_name_dict = dict(zip(meme_motif_df['motif_id'], meme_motif_df['dedup_protein_name']))
    meme_motif_df = meme_motif_df[meme_motif_df['pval'] < pco]
    meme_motif_df = meme_motif_df.value_counts(['peak', 'dedup_protein_name']).unstack()
    meme_motif_df.index = meme_motif_df.index.str.lower().str.replace(":", "_").str.replace("-", "_")
    meme_motif_df = meme_motif_df.fillna(0)
    meme_motif_df = mod_index(meme_motif_df)
    
    vs = []
    missing_inds = atac_peaks_index
    for i in missing_inds:
        v = pd.Series(index=meme_motif_df.columns).fillna(0)
        vs.append(v)
    
    tmpdf = pd.concat(vs, axis=1).T
    tmpdf.index = missing_inds
    meme_motif_df = pd.concat([meme_motif_df, tmpdf], axis=0)
    meme_motif_df = meme_motif_df.loc[atac_peaks_index.index]
    return meme_motif_df, motif_id_to_name_dict


def make_metadomain_hub_freq_df(self, inter_and_intra_connections_tcon, inter_and_intra_connections_treg):
    metadomain_mat_dict = {
            'Treg' : (inter_and_intra_connections_tcon == 0) & (inter_and_intra_connections_treg > 0),
            'Both' : (inter_and_intra_connections_tcon > 0) & (inter_and_intra_connections_treg > 0),
            'Tcon' : (inter_and_intra_connections_tcon > 0) & (inter_and_intra_connections_treg == 0),
            'Neither' : (inter_and_intra_connections_tcon == 0) & (inter_and_intra_connections_treg == 0),
    
    }
    
    metadomain_hub_freq_df = pd.DataFrame()
    for key, metadomain_mat in metadomain_mat_dict.items():
        data = []
        us = [0, 4, 18]
        for u in us:
            inds = self.goodinds[self.merged_clustdict['all']==u]
            data.append(np.sum(metadomain_mat[inds, :][:, inds]))
        series_data = pd.Series(data, index=us)
        metadomain_hub_freq_df[key] = series_data
    metadomain_hub_freq_df.index = [columns_to_names.get(x, "Other") for x in metadomain_hub_freq_df.index]
    metadomain_hub_freq_df = (metadomain_hub_freq_df.T / metadomain_hub_freq_df.sum(axis=1)).T
    metadomain_hub_freq_df = metadomain_hub_freq_df.drop("Neither", axis=1)    
    return metadomain_hub_freq_df

def load_bulkrna(my_tss_df, geneLengths):
    active_gene_lfcs = pbt.BedTool('./gene_expression/Treg_actv_vs_Tcon_actv_thresh=0.25.csv.narrowPeak').to_dataframe().set_index('name')
    resting_gene_lfcs = pbt.BedTool('./gene_expression/Treg_rest_vs_Tcon_rest_thresh=0.25.csv.narrowPeak').to_dataframe().set_index('name')

    active_gene_lfcs = active_gene_lfcs[active_gene_lfcs['itemRgb'] > 4]
    resting_gene_lfcs = resting_gene_lfcs[resting_gene_lfcs['itemRgb'] > 4]

    active_gene_lfcs = active_gene_lfcs[active_gene_lfcs.index.isin(my_tss_df['gene_name']) & 
                                        active_gene_lfcs.index.isin(geneLengths.index)]
    resting_gene_lfcs = resting_gene_lfcs[resting_gene_lfcs.index.isin(my_tss_df['gene_name']) & 
                                        resting_gene_lfcs.index.isin(geneLengths.index)]

    active_gene_lfcs = active_gene_lfcs[~active_gene_lfcs.index.duplicated()]
    resting_gene_lfcs = resting_gene_lfcs[~resting_gene_lfcs.index.duplicated()]

    active_gene_lfcs['rpkm'] = active_gene_lfcs['itemRgb'] / geneLengths.loc[active_gene_lfcs.index]
    resting_gene_lfcs['rpkm'] = resting_gene_lfcs['itemRgb'] / geneLengths.loc[resting_gene_lfcs.index]

    gene_dict = {
        'Active' : active_gene_lfcs,
        'Resting' : resting_gene_lfcs,
    }    
    return gene_dict