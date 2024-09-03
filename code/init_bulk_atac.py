import glob
import pybedtools as pbt
# from init_bulk_ChIP import all_peak_dict
import sys
sys.path.append('../../code/')
from aux_functions import get_col, remove_chr_bedtool
import pandas as pd
atac_dict = {}
for file in glob.glob('./atac/processed/*thresh=0*.csv'):
    cond1, cond2 = file.split("/")[-1].split("_thresh")[0].split("_vs_")
    thresh = file.split(".csv")[0].split("thresh=")[1]
    atac_dict[f'{cond1}_vs_{cond2}_thresh={thresh}'] = pbt.BedTool(file)


def get_index_from_bedtool_df(df):
    ind = df['chrom'].astype(str) + "_" + df['start'].astype(str) + "_" + df['end'].astype(str) 
    return ind


atac_peak_bedtool = atac_dict['Treg_actv_vs_Tcon_actv_thresh=0']
atac_peak_df = atac_peak_bedtool.to_dataframe()
ind = get_index_from_bedtool_df(atac_peak_df)
atac_peak_df.index = ind
