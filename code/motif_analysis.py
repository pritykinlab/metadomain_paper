import scipy
import scipy.ndimage
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import pybedtools as pbt


def nochr_make_ctcf_peaks(anchors, name, shift=0, outdir='./meme/for_motifs/', fapath='./annotations/mm10.fa', genome='mm10'):
    fasta = pbt.BedTool(fapath)
    coords = anchors.shift(genome=genome, s=shift)
    len(coords)
    seq = open(coords.sequence(fi=fasta).seqfn).read()
    with open(f'{outdir}/{name}', 'w+') as f:
        f.write(seq.upper())
    return seq
	

import time

from sklearn.preprocessing import LabelEncoder

import time
def grab_fimo_output_and_make_motif_df(placedict, namedict, shiftdict = {}):
    all_possible_places = []
    all_possible_tuples = []
    for x in placedict.values():
        for y in x:
            tup = tuple(list(y)[:3])
            place = '_'.join(list(y)[:3])
            all_possible_places.append(place)
            all_possible_tuples.append(tup)
    boundary_enc = LabelEncoder()
    boundary_enc.fit(all_possible_places)
    bedtool_dict = {}
    for name in placedict:
        t1 = time.time()
        print(f'./meme/fimo_output/{name}/fimo.tsv', f'./peaks/motifs/{name}.bed')
        motif_bedtool = pbt.BedTool(f'./peaks/motifs/{name}.bed')
        bedtool_dict[name] = motif_bedtool
        t1f = time.time()
    print("Done with fimo to bed", t1f-t1)

    all_motif_names = set()
    for name, motif_bed in bedtool_dict.items():
        newnames = set(get_col(motif_bed, 3))
        all_motif_names = all_motif_names.union(newnames)

    motif_enc = LabelEncoder()
    motif_enc.fit(list(all_motif_names))

    dfdict = {}
    for name, bed in bedtool_dict.items():
        motif_df = make_motif_df(pbt.BedTool(all_possible_tuples), all_motif_names, motif_enc, boundary_enc,
                             file=bed.fn, shift=shiftdict.get(name, 0))
        motif_df.columns = [namedict[x] for x in motif_df.columns]
        dfdict[name] = motif_df
    return dfdict


def fimo_to_bed(file, out):
    with open(file) as f:
        print_lines = []
        for i, line in enumerate(f):
            t1 = time.time()
            if i==0:
                continue
            elif line=='\n':
                break
            elif line[0] == '#':
                continue
            else:
                vals = line.split('\t')
                motif_name = vals[0]
                protein_name = vals[1]
                chrom, _ = vals[2].split(':')
                peak = vals[2]
                loop_start, loop_end = _.split('-')
                loop_start, loop_end = map(int, [loop_start, loop_end])
                motif_start, motif_end = map(int, [vals[3], vals[4]])
                sequence = vals[-1]
                start = loop_start+motif_start
                end = loop_start+motif_end
                strand = vals[5]
                score = float(vals[6])
                pval = float(vals[7])
                qval = float(vals[8])  
                print_lines.append([str(chrom), str(start), str(end), str(peak), str(motif_name), str(protein_name), str(strand), str(score), str(pval), str(qval), str(sequence)])
            t1f = time.time()
    with open(out, 'w+') as f:
        for line in print_lines:
            f.write("\t".join(line) + "\n")

from collections import Counter
def fimo_to_df(file):
    with open(file) as f:
        rows = []
        for i, line in enumerate(f):
            t1 = time.time()
            if i==0:
                continue
            elif line=='\n':
                break
            else:
                vals = line.split('\t')
                motif_name = vals[0]
                protein_name = vals[1]
                chrom, _ = vals[2].split(':')
                peak = vals[2]
                loop_start, loop_end = _.split('-')
                loop_start, loop_end = map(int, [loop_start, loop_end])
                motif_start, motif_end = map(int, [vals[3], vals[4]])
                start = loop_start+motif_start
                end = loop_start+motif_end
                strand = vals[5]
                score = float(vals[6])
                pval = float(vals[7])
                qval = float(vals[8])  
                rows.append([str(chrom), str(start), str(end), str(peak), str(motif_name), str(protein_name), str(strand), float(score), float(pval), float(qval), ])
    fimo_motif_df = pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'peak', 'motif_name', 'protein_name', 'strand', 'score', 'pval', 'qval'])
    l = fimo_motif_df.value_counts(['protein_name', 'motif_name']).reset_index()

    c = Counter()
    name_to_protein = {}
    for _, row in l.iterrows():
        protein, name = row.protein_name, row.motif_name
        if c[protein] > 0 :
            name_to_protein[name] = f'{protein}.{c[protein]}'
        else:
            name_to_protein[name] = f'{protein}'
        c[protein] += 1
    fimo_motif_df['protein_name_dedup'] = fimo_motif_df['motif_name'].apply(name_to_protein.get)
    return fimo_motif_df

arr = np.asarray
def get_col(bed, col):
    return arr(list(zip(*bed))[col])

def get_first_n_columns(bedtool, n=3):
    new_bedtool = pbt.BedTool.from_dataframe(bedtool.to_dataframe().iloc[:, :n])
    return new_bedtool

def make_motif_df(boundaries_as_bedtool, all_motif_names, motif_enc, boundary_enc, file='./peaks/motifs/all_anchor_motifs.bed', shift=0):
    n_motifs = len(all_motif_names); n_ancs = len(boundaries_as_bedtool)
    anc_motif_mat = np.zeros((n_ancs, n_motifs))
    anclist = []; motiflist = []
    sub_bedtool = get_first_n_columns(boundaries_as_bedtool, n=3)
    for i in sub_bedtool.intersect(pbt.BedTool(file).shift(s=shift, g=chromsizepath), wo=True):
        anc = i[:3]
        motif = i[6]
        anclist.append('_'.join(anc))
        motiflist.append(motif)

    i1s = boundary_enc.transform(anclist)
    i2s = motif_enc.transform(motiflist)
    for i1, i2 in zip(i1s, i2s):
        anc_motif_mat[i1, i2] += 1

    motif_df = pd.DataFrame(anc_motif_mat)
    motif_df.columns = motif_enc.inverse_transform(motif_df.columns.values)
    motif_df.index = boundary_enc.inverse_transform(motif_df.index)
    return motif_df

import statsmodels
import statsmodels.stats
import statsmodels.stats.multitest
def get_stats_from_motif_dfs(df1, df2):
    stat_df = []
    for col in df1.columns:
        v1, v2 = df1[col], df2[col]
        stat, pval = scipy.stats.ranksums(v1, v2)
        delta = np.log2(np.mean(v1) / np.mean(v2))
        stat_df.append([col, stat, pval, delta])
    df = pd.DataFrame(stat_df, columns=['Motif', 'Z', 'pval', 'stat'])
    df['padj'] = statsmodels.stats.multitest.fdrcorrection(df['pval'])[1]
    return df



def get_motif_enrichments_delta(motif_df, df1, df2):
    rs = []
    for c, col in enumerate(motif_df):
        motif_idx = motif_df[col] > 0

        x = scipy.stats.zscore(df1.loc[:, 'log2FoldChange']).loc[motif_idx]
        y = scipy.stats.zscore(df2.loc[:, 'log2FoldChange']).loc[motif_idx]
        _, p = scipy.stats.ranksums(x, y)
        delta = (x.mean() - y.mean())
        rs.append([col, p, delta])    
    return pd.DataFrame(rs, columns = ['factor', 'p', 'delta'])






