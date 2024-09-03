import cooler 
import numpy as np
import pandas as pd
import pybedtools as pbt
import matplotlib.pyplot as plt

arr = np.asarray

def make_rna_cdf_plots(tdict, peak_genes, key_to_name):
    fig, ax = plt.subplots(figsize=(4, 3))
    for cond, t in tdict.items():
        g0 = peak_genes.intersect(t, u=True)
        pvals = get_col(g0, 4, float)
        vs = []
        pcos = np.log(np.geomspace(1e-20, .05, 1000))
        for pco in pcos:
            vs.append((np.log(pvals) < pco).mean())
        ax.plot(pcos, vs, label=cond)
    ax.legend()

    vdict = {}
    fig, ax = plt.subplots(1, figsize=(5, 4))
    for cond, t in tdict.items():
        g0 = peak_genes.intersect(t, u=True)
        lfcs = (get_col(g0, 6, float))
        pvals = np.abs(get_col(g0, 4, float))
        vs = []
        pcos = (np.linspace(np.min(lfcs)//1-1, np.max(lfcs)//1 + 1, 1000))
        for pco in pcos:
            vs.append((lfcs < pco).mean())
        if 'all' not in cond:
            ax.plot(pcos, vs, label=f'{key_to_name[cond]} > 0')
        else:
            ax.plot(pcos, vs, label=f'All')        
        vdict[cond] = lfcs
    ax.set_xlabel("Tcon ÷ Treg RNA-seq")
    ax.set_ylabel("Proportion")

    ax.legend(bbox_to_anchor=(1, 1))


def get_col(bed, col, dtype=None):
    if dtype is None:
        return arr(list(zip(*bed))[col])
    else:
        return arr(list(zip(*bed))[col]).astype(dtype)

import scipy
def make_gene_comparison_plot(w1, w2, genes, labels, title, xlabel, cdict= {}, xmin=-6, xmax=6, lw=10):
    overlap = []
    overlap_bt = w1.intersect(w2, u=True)
    for i in genes.intersect(overlap_bt, u=True):
        overlap.append(float(i[6]))
    w1_alone = []
    w1_alone_bt = w1.subtract(w2, A=True)
    for i in genes.intersect(w1_alone_bt, u=True):
        w1_alone.append(float(i[6]))

    w2_alone = []
    w2_alone_bt = w2.subtract(w1, A=True)
    for i in genes.intersect(w2_alone_bt, u=True):
        w2_alone.append(float(i[6]))
    all_fcs = arr(list(zip(*genes))[6]).astype(float)

    fig, ax = cdf([overlap, w1_alone, w2_alone, all_fcs], labels, xmin=xmin, xmax=xmax, cdict=cdict, linewidth=lw)
    ax.legend(bbox_to_anchor = (1, 1))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    return overlap_bt, w1_alone_bt, w2_alone_bt, fig, ax






def make_closest_gene_comparison_plot(w1, w2, peak_genes, labels, title, xlabel, xmin=-6, xmax=6, cdict={}):
    overlap = []
    overlap_bt = w1.intersect(w2, u=True)
    test = overlap_bt.closest(peak_genes.sort(), d=True).filter(lambda x: int(x[-1]) >= 0).saveas()
    new_genes = []
    for i in test:
        new_genes.append(tuple(i[6:]))
    new_genes = pbt.BedTool(list(set(new_genes)))
    for i in new_genes:
        overlap.append(float(i[6]))

    w1_alone = []
    w1_alone_bt = w1.subtract(w2, A=True)
    test = w1_alone_bt.closest(peak_genes.sort(), d=True).filter(lambda x: int(x[-1]) >= 0).saveas()
    new_genes = []
    for i in test:
        new_genes.append(tuple(i[6:]))
    new_genes = pbt.BedTool(list(set(new_genes)))
    for i in new_genes:
        w1_alone.append(float(i[6]))

    w2_alone = []
    w2_alone_bt = w2.subtract(w1, A=True)
    test = w2_alone_bt.closest(peak_genes.sort(), d=True).filter(lambda x: int(x[-1]) >= 0).saveas()
    new_genes = []
    for i in test:
        new_genes.append(tuple(i[6:]))
    new_genes = pbt.BedTool(list(set(new_genes)))
    for i in new_genes:
        w2_alone.append(float(i[6]))
    all_fcs = arr(list(zip(*peak_genes))[6]).astype(float)

    fig, ax = cdf([overlap, w1_alone, w2_alone, all_fcs], labels, xmin=xmin, xmax=xmax, cdict=cdict)
    ax.legend(bbox_to_anchor = (1, 1))
    ax.set_title(title)
    ax.set_xlabel(xlabel)


import itertools
def cdf(fcs, names, cdict = {}, xlabel="", folder="tmp", filename="tmp", xmin=None, xmax=None, **kwargs):

    lens ={}
    if xmin is None:
        xmin = np.min(fcs[0])
    if xmax is None:
        xmax = np.max(fcs[0])
    xs = np.linspace(xmin, xmax, 1000)
    fc_dict = {}
    for i, fc in enumerate(fcs):
        name = names[i]
        fc_dict[name] = []
        for j, x in enumerate(xs):
            fc_dict[name].append((np.asarray(fc) < x).mean())
        lens[name] = len(fc)
    fig, axs = plt.subplots(figsize=(5, 3))
    for name in fc_dict:
        axs.plot(xs, fc_dict[name], label=f'{name}, {lens[name]}', color=cdict.get(name), **kwargs)
    axs.legend()

    axs.set_xlabel(xlabel)
    axs.set_xlim([xmin, xmax])
    plt.tight_layout()

    dict_for_ks = dict(zip(names, fcs))
    for pair in itertools.combinations(list(dict_for_ks.keys()), 2):
        name1, name2 = pair
        ks_test = scipy.stats.ks_2samp(dict_for_ks[name1], dict_for_ks[name2])
        pval = ks_test[1]
        print(pval)

    return fig, axs







