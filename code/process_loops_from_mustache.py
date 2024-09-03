import numpy as np
import pandas as pd
import pybedtools as pbt
from aux_functions import tuple_to_grange, remove_chr_bedtool

def loops_from_df(df, pco=.05):
    return pbt.BedTool.from_dataframe(df)


def remove_chr_bedtool_loops(bedtool):
    new_bedtool_list = []
    for i in bedtool:
        z = list(i) 
        chrom = z[0]
        if 'chr' == chrom[:3]:
            newchrom = chrom[3:]
        else:
            newchrom = chrom
        z[0] = newchrom

        chrom = z[3]
        if 'chr' == chrom[:3]:
            newchrom = chrom[3:]
        else:
            newchrom = chrom
        z[3] = newchrom
        new_bedtool_list.append(z)
    new_bedtool = pbt.BedTool(new_bedtool_list)
    return new_bedtool

def filter_raw_mustache_loops(mustache_loop_df):
    mustache_loop_df = mustache_loop_df.copy()
    pbt_loops = remove_chr_bedtool_loops(loops_from_df(mustache_loop_df))
    loops = {}
    for loop in pbt_loops:
        l1 = tuple(loop[:3])
        l2 = tuple(loop[3:6])
        loops[(l1, l2)] = 1    

    ancs = {}
    for l in pbt_loops:
        l1, l2 = l[:3], l[3:6]
        ancs[(l1[0], int(l1[1]), int(l1[2]))]=1
        ancs[(l2[0], int(l2[1]), int(l2[2]))]=1

    res=5000
    sorted_ancs = sorted(list(ancs))
    anc_mapping = {}
    merge_ancs = []
    c=1
    for i, anc in enumerate(sorted_ancs[:-1]):
        next_anc = sorted_ancs[i+1]
        if (anc[0] == next_anc[0]) and (np.abs(anc[1] - next_anc[1]) < 5000*3):
            merge_ancs.append(anc)
            c+=1
        else:
            if c > 1:
                merge_ancs.append(anc)
                if c%2 != 0:
                    start = np.median(list(zip(*merge_ancs))[1])
                else:
                    starts = list(zip(*merge_ancs))[1]
                    midpoint = np.mean(starts)//5000*5000
                    start = int(midpoint)
                for a in merge_ancs:
                    v = a[0], int(start), int(start+res)
                    assert v[2] - v[1] > 4000
                    anc_mapping[a] = (v)
                merge_ancs = []
            else:
                new_anc = list(anc)
                new_anc[2] = new_anc[1] + 5000
                anc_mapping[anc] = new_anc
                if i==len(sorted_ancs)-2:
                    anc_mapping[next_anc] = next_anc
            c=1
    newloopset = remap_ancs(pbt_loops, anc_mapping)
    modified_df = pbt_loops.to_dataframe()
    modified_df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'fdr', 'size']
    return newloopset, anc_mapping, modified_df

def remap_ancs(pbt_loops, anc_mapping):
    newloopset = set()
    for loop in pbt_loops:
        l1, l2 = loop[:3], loop[3:6]
        l1 = (l1[0], int(l1[1]), int(l1[2]))
        l2 = (l2[0], int(l2[1]), int(l2[2]))
        try:
            l1_new = tuple(anc_mapping.get(l1))
            l2_new = tuple(anc_mapping.get(l2))
        except Exception:
            print(l1, l2)
            raise Exception
        newloopset.add(l1_new + l2_new)
    return pbt.BedTool(list(newloopset))


def loopset_to_granges(loopset):
    grange_list = []
    for loop in loopset:
        l1, l2 = loop[:3], loop[3:6]
        loop_grange = tuple_to_grange(*l1) + "|" + tuple_to_grange(*l2)
        grange_list.append(loop_grange)
    return grange_list


def add_loops_to_reference(reference_loops, loops_to_add, res=5000, d = 2):
    remapped_loops_to_add = loops_to_add.copy()
    
    reference_loop_df = reference_loops.copy()
    loops_to_add_df = loops_to_add.copy()
    
    same_chr = np.equal.outer(reference_loop_df['BIN1_CHR'].values, loops_to_add_df['BIN1_CHR'].values)
    dist_in_x = np.subtract.outer(reference_loop_df['BIN1_START'].values, loops_to_add_df['BIN1_START'].values)
    dist_in_y = np.subtract.outer(reference_loop_df['BIN2_START'].values, loops_to_add_df['BIN2_START'].values)

    close_in_x = np.abs(dist_in_x) <= res*d
    close_in_y = np.abs(dist_in_y) <= res*d
    has_reference_loop_matrix = (close_in_x & close_in_y & same_chr)
    have_reference_loop = has_reference_loop_matrix.any(axis=0)
    # # # Remap the loops which have a reference loop
    map_to_reference = {}
    for loop_ind in np.where(have_reference_loop)[0]:
        # print(loop_ind)
        reference_loop_ind = np.where(has_reference_loop_matrix[:, loop_ind])[0][0]
        map_to_reference[loop_ind] = reference_loop_ind
    for loop_ind, reference_loop_ind in map_to_reference.items():
        remapped_loops_to_add.iloc[loop_ind] = reference_loop_df.iloc[reference_loop_ind].copy()
        
    assert have_reference_loop.shape[0] == len(loops_to_add_df)
    loops_to_add_df = loops_to_add_df[~have_reference_loop]
    reference_loop_df = pd.concat([reference_loop_df, loops_to_add_df], axis=0)
    return reference_loop_df, remapped_loops_to_add

