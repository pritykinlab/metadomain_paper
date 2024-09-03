import cooler
import numpy as np

def initialize_helper_vars(coolfile,  parsed_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
                     '12', '13', '14', '15', '16', '17', '18', '19', 'X']):
    chromsizes = {}
    for chrom, size in coolfile.chromsizes.items():
        chromsizes[chrom] = size

    inds_to_region = {}
    region_to_inds = {}
    for chrom in parsed_chroms:
        df = coolfile.bins()[:]
        subdf = df[df['chrom'] == chrom]
        regions = list(zip(subdf['chrom'].values, map(int, subdf['start'].values), map(int, subdf['end'].values)))
        inds = np.arange(len(regions))
        inds_to_region[chrom] = regions
        region_to_inds[chrom] = dict(zip(regions, inds))

    all_region_to_ind = {}
    tot = 0
    for chrom in parsed_chroms:
        df = coolfile.bins()[:]
        subdf = df[df['chrom'] == chrom]
        regions = list(zip(subdf['chrom'].values, map(int, subdf['start'].values), map(int, subdf['end'].values)))
        inds = np.arange(len(regions)) + tot
        for r, i in zip(regions, inds):
            all_region_to_ind[r] = i
        tot += len(regions)
    all_ind_to_region = list(all_region_to_ind.keys()) 
    chrom_to_start = {}
    chrom_to_end = {}
    c = 0 
    for _, chrom in enumerate(parsed_chroms):
        if ('M' in chrom) or ('Y' in chrom):
            continue
        vs = coolfile.matrix(balance=False).fetch(chrom)
        n = len(vs)
        chrom_to_start[chrom] = c
        c += n
        chrom_to_end[chrom] = c   
    return chromsizes, parsed_chroms, region_to_inds, all_region_to_ind, inds_to_region, all_ind_to_region, chrom_to_start, chrom_to_end

def make_int(l):
   return l[0], int(l[1]), int(l[2])

import pybedtools as pbt
from aux_functions import add_chr_to_bedtool
def initialize_genes(all_ind_to_region, all_region_to_ind, path_to_genes='./peaks/RNA_coverage.narrowPeak',
                     filetype=None):
    if filetype == 'GTF':
        gene_gtf = pbt.BedTool(path_to_genes)
        for i in gene_gtf.intersect(all_ind_to_region, wo=True):
            gene = i[3]
            if "Rik" in gene:
                continue
            tad = tuple(i[-4:-1])
            ind = all_region_to_ind[make_int(tad)]
            gene_to_ind.setdefault(gene, [])
            gene_to_ind[gene].append(ind)
            pval = float(i[4])
            fc = float(i[6])
            ind_to_gene.setdefault(ind, [])
            ind_to_gene[ind].append(gene)
            ind_to_fcs.setdefault(ind, [])
            ind_to_fcs[ind].append(fc)
            ind_to_pvals.setdefault(ind, [])
            ind_to_pvals[ind].append(pval)

    else:            
        if 'chr' in all_ind_to_region[0][0]:
            peak_genes = add_chr_to_bedtool(pbt.BedTool(path_to_genes))
        else:
            peak_genes = pbt.BedTool(path_to_genes)

        gene_to_ind = {}
        ind_to_gene = {}
        ind_to_fcs = {}
        ind_to_pvals = {}
        for i in peak_genes.intersect(all_ind_to_region, wo=True):
            gene = i[3]
            if "Rik" in gene:
                continue
            tad = tuple(i[-4:-1])
            ind = all_region_to_ind[make_int(tad)]
            gene_to_ind.setdefault(gene, [])
            gene_to_ind[gene].append(ind)
            pval = float(i[4])
            fc = float(i[6])
            ind_to_gene.setdefault(ind, [])
            ind_to_gene[ind].append(gene)
            ind_to_fcs.setdefault(ind, [])
            ind_to_fcs[ind].append(fc)
            ind_to_pvals.setdefault(ind, [])
            ind_to_pvals[ind].append(pval)
        return gene_to_ind, ind_to_gene

