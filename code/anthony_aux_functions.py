from aux_functions import index_to_bedtool, bedtool_to_index
def overlap_two_separate_dfs(df1, df2):
    bed1, bed2 = index_to_bedtool(df1.index), index_to_bedtool(df2.index)
    l1s = []
    l2s = []
    for i in bed1.intersect(bed2, wo=True):
        l1 = i[:3]
        l2 = i[3:6]
        l1s.append(l1)
        l2s.append(l2)
    return bedtool_to_index(l1s), bedtool_to_index(l2s)

def df_bedtool_overlaps(df, bedtool):
    overlaps = pd.Series(get_col(index_to_bedtool(df.index).intersect(bedtool, c=True), -1).astype(int),
                         index = df.index)
    return overlaps



# xiao_pref=/Genomics/levineshare/People/Xiao/trl_mcools
# files=(halo_trl_merge.mcool delta_poz_merge.mcool)
# for file in ${files[@]}
# do 
#     cooltools expected-cis ${xiao_pref}/${file}::/resolutions/400 -p 16 -o  expected/expected_xiao_${file}.tsv
# done
