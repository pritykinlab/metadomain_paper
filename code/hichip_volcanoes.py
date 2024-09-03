import numpy as np
import scipy
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
import matplotlib.pyplot as plt
import bbi
import pybedtools as pbt
import pickle
from adjustText import adjust_text


def make_text_for_volcano(x, y, indsoi, xinds, yinds, ax, upcolor='orange', downcolor='orange', scatter=True, ind_to_gene = {}, expand_text=(.2, .2), expand_points=(.2, .2), 
				s=80, fontsize=12):
    texts = []
    for i in indsoi:
        xind = xinds[i]; yind = yinds[i]
        name1 = get_name(xind, ind_to_gene)
        name2 = get_name(yind, ind_to_gene)
        name = f"{name1} to {name2}"
        if 'None' in name:
                continue
        if x[i] < 0:
            t = ax.text(x[i], y[i], name, ha = 'left', va='bottom', fontsize=fontsize, zorder=10, bbox=dict(facecolor='none', edgecolor='none', pad=0))
        else:
            t = ax.text(x[i], y[i], name, ha = 'right', va='bottom', fontsize=fontsize, zorder=10, bbox=dict(facecolor='none', edgecolor='none', pad=0))
        if scatter:
            if x[i] > 0:
                color=upcolor
            else:
                color=downcolor
            ax.scatter(x[i], y[i], s=s, color=color, edgecolor='black', linewidths=1, zorder=12)

        texts.append(t)
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1,  ),
        expand_text=(1.51, 1.55), expand_points=(1.51, 1.55),
               force_text = expand_text, force_points = expand_points, zorder=1)    

def get_names_from_loop(loop, anc_to_name):
    l1, l2 = loop[:3], loop[3:6]
    l1, l2 = tuple(l1), tuple(l2)
    names1, names2 = anc_to_name.get(l1, ['None']), anc_to_name.get(l2, ['None'])
    return names1, names2 

def make_text_for_loop_volcano(x, y, indsoi, loopsoi, ax, upcolor='orange', downcolor='orange', 
                scatter=True, anc_to_name = {}, expand_text=(.2, .2), expand_points=(.2, .2), s=80, fontsize=12):
    texts = []
    for i in indsoi:
        i = int(i)
        loop = loopsoi[i]
        names1, names2 = get_names_from_loop(loop, anc_to_name)
        name = f"{names1[0]} to {names2[0]}"
        if ('None' in names1) and ('None' in names2):
            continue
        elif ('None' in names1):
            name = f"{names2[0]}"
        elif ('None' in names2):
            name = f"{names1[0]}"

        if x[i] < 0:
            t = ax.text(x[i], y[i], name, ha = 'left', va='bottom', fontsize=fontsize, zorder=10, bbox=dict(facecolor='none', edgecolor='none', pad=0))
        else:
            t = ax.text(x[i], y[i], name, ha = 'right', va='bottom', fontsize=fontsize, zorder=10, bbox=dict(facecolor='none', edgecolor='none', pad=0))
        if scatter:
            if y[i] > 0:
                color=upcolor
            else:
                color=downcolor
            ax.scatter(x[i], y[i], s=s, color=color, edgecolor='black', linewidths=1, zorder=12)

        texts.append(t)
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1,  ),
        expand_text=(1.51, 1.55), expand_points=(1.51, 1.55),
               force_text = expand_text, force_points = expand_points, zorder=1)    


def make_text_for_1d_volcano(x, y, indsoi, inds, ax, upcolor='orange', downcolor='orange', scatter=True, ind_to_gene = {}, expand_text=(.2, .2), expand_points=(.2, .2),
                                s=80):
    texts = []
    for i in indsoi:
        ind = inds[i];
        name1 = get_name(ind, ind_to_gene)
        name = f"{name1}"
        if 'None' in name:
                continue
        if x[i] < 0:
            t = ax.text(x[i], y[i], name, ha = 'left', va='bottom', fontsize=10, zorder=10, bbox=dict(facecolor='none', edgecolor='none', pad=0))
        else:
            t = ax.text(x[i], y[i], name, ha = 'right', va='bottom', fontsize=10, zorder=10, bbox=dict(facecolor='none', edgecolor='none', pad=0))
        if scatter:
            if x[i] > 0:
                color=upcolor
            else:
                color=downcolor
            ax.scatter(x[i], y[i], s=s, color=color, edgecolor='black', linewidths=1, zorder=2)
        texts.append(t)
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1,  ),
        expand_text=(1.51, 1.55), expand_points=(1.51, 1.55),
               force_text = expand_text, force_points = expand_points, zorder=0)








def unzreg(l):
    chrom, s, e = list(zip(*l))
    return chrom, list(arr(s).astype(int)), list(arr(e).astype(int))


arr = np.asarray
from copy import deepcopy
from matplotlib.markers import MarkerStyle
from mpl_toolkits import axes_grid1
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def make_hichip_volcano(sub_df, cond, hic_cond, ind_to_gene, all_connections, all_ind_to_region, foxp3_bw, cond2name, FDR_CO = .1, maxy=2, **kwargs):
    indsoi = np.where(np.sum(all_connections, axis=1)>0)[0]
    chrom, s, e = unzreg(all_ind_to_region[x] for x in indsoi)
    all_vs = foxp3_bw.stackup(chrom, s, e, bins=1)
    avg_vs = np.nanmean(all_vs)
    x, y = sub_df['log2FoldChange'].values, sub_df['newpadj'].values
    x, y = deepcopy(x), deepcopy(y)
    y = -np.log10(y)
    xinds, yinds = zip(*[x.split("_")  for x in sub_df.index])
    xinds, yinds = list(map(int, xinds)), list(map(int, yinds))

    chrom, s, e = unzreg([all_ind_to_region[x] for x in xinds])
    vs = foxp3_bw.stackup(chrom, s, e, bins=1)
    chrom, s, e = unzreg([all_ind_to_region[x] for x in yinds])
    vs2 = foxp3_bw.stackup(chrom, s, e, bins=1)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=60, c=vs, edgecolors='gray', zorder=1, vmin=1, vmax=3,
               cmap=cm.bwr, marker=MarkerStyle('o', fillstyle='left'), label='FoxP3 signal, left bin'
              )

    pos = ax.scatter(x, y, s=60, c=vs2, edgecolors='gray', zorder=1, cmap=cm.bwr, vmin=1, vmax=2,
               marker=MarkerStyle('o', fillstyle='right'), label='FoxP3 signal, right bin', 
              )

    ax.legend(bbox_to_anchor=(1.15, 1), loc='upper left')
    indsoi = (np.argsort(-y)[:10])
    make_text_for_volcano(x, y, indsoi, xinds, yinds, ax=ax, scatter=False, ind_to_gene=ind_to_gene, **kwargs)
    
    max_x_val = np.max(np.abs(y))*1.1
    max_x_val = max(1.5, max_x_val)
    ax.set_xlim([-max_x_val, max_x_val])
    ax.set_ylim([-.1, maxy])
    ax.hlines(-np.log10(FDR_CO), -max_x_val, max_x_val, color='gray', linestyle='--',)
    ax.text(max_x_val, -np.log10(FDR_CO), s=f'FDR={FDR_CO}', ha = 'right', va='bottom', fontsize=12)
    
    ax.set_ylabel("-log10(FDR)")
    ax.set_xlabel(f"LFC ({cond2name.get(cond)})")
    ax.set_title(f"Differential HiChIP at \nfocal contacts where {hic_cond}")

    newax = inset_axes(ax, width = '2%', height = '100%', 
                        loc = 'lower left', bbox_to_anchor=(1.1, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
    plt.colorbar(pos, cax=newax)
    newax.set_yticks([1, avg_vs, 2])
    newax.set_yticklabels([1, "Average in bins w/ focal contacts", 2, ])
    df = pd.DataFrame()
    df['vL'] = vs.flatten()
    df['vR'] = vs2.flatten()
    df['x'] = x
    df['y'] = y
    df['xinds'] = xinds
    df['yinds'] = yinds
    return df 







