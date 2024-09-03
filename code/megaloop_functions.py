from aux_functions import *

def get_megaloops_by_distance(megaloop_mat):
    x, y = np.where(megaloop_mat)
    mega_df = pd.DataFrame([x, y]).T
    mega_df.columns = ['row', 'col']
    mega_df['dist'] = (mega_df['col'] - mega_df['row']).abs()
    mega_df['score'] = megaloop_mat[x, y]
    mega_df = mega_df[mega_df['row'] < mega_df['col']]
    mega_df['id'] = mega_df['row'].astype(str) + "_" +  mega_df['col'].astype(str)
    return mega_df

def subset_ind_mega_df_by_column(mega_df, clow, chigh, col='dist'):
    mega_df = mega_df[(mega_df[col] > clow) & (mega_df[col] < chigh)]
    return mega_df

def ind_mega_df_to_set(ind_mega_df):
    vals = set(ind_mega_df['row'].astype(str) + '_' + ind_mega_df['col'].astype(str))
    return vals

def ind_mega_df_to_set_by_col(ind_mega_df, clow, chigh, col='dist'):
    ind_mega_df = subset_ind_mega_df_by_column(ind_mega_df, clow, chigh, col=col)
    vals = ind_mega_df_to_set(ind_mega_df)
    return vals
