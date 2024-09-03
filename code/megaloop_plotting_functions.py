import seaborn as sns
import scipy
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
from copy import deepcopy
from plotting_functions import add_bigwig_to_axis
from plotting_functions import add_GTF_to_axis
from plotting_functions import add_GTF_to_L_axis
from aux_functions import extend_l
def make_scoping_plot_for_fig1(coolfile, row, col, get_func, name_of_func, vert_bbdict = {},
                                resolution=5000, d = 110000, 
                                d_RR=None, 
                                d_LL=None, 
                                d_RL=None,
                                d_LR=None,
                                useSigma=False, ignore_set=[],
                                vmin=0, vmax=1, balance=True, bwcolor='black', 
                                ylimdict=None, all_ind_to_region=None, 
                                filesoi=None, 
                                ax = None,
                                dpi=100, ignore_method='ignore',
                                **kwargs):
    if d_RR is None:
        d_RR = d
    if d_LL is None:
        d_LL = d
    if d_RL is None:
        d_RL = d
    if d_LR is None:
        d_LR = d
    if type(row) == int:
        chrom1, s1, e1 = all_ind_to_region[row]
        place1  = (chrom1, s1-d_LL, e1+d_LR)

        chrom2, s2, e2 = all_ind_to_region[col]
        place2  = (chrom2, s2-d_RL, e2+d_RR)
    else:
        chrom1, s1, e1 = row
        chrom2, s2, e2 = col
        place1, place2 = row, col
        place1, place2 = extend_l(row, d), extend_l(col, d)

    m, _, _ = get_func(coolfile, place1, place2, filesoi=filesoi, resolution=resolution, balance=balance)
    fullm = m
    print('SHAPE', fullm.shape)
    
    sns.set_style('ticks')
    if useSigma:
        m2 = m.copy()
        m2[np.isnan(m2)] = 0
        filtered_mat = scipy.ndimage.gaussian_filter(m2, sigma=.75)
    else:
        filtered_mat = fullm
        print(np.max(filtered_mat))
    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 3), dpi=dpi)
    else:
        fig = None
    if (place2[1] <= place1[2]) and (chrom1 == chrom2) and (place2[2] >= place1[2]):
        ax.matshow(filtered_mat, cmap=cm.gist_heat_r, norm=colors.LogNorm(vmin=vmin, vmax=vmax), extent=[place2[1], place2[2], place1[2], place1[1], ]
                  , rasterized=True)
    else:
        ax.matshow(filtered_mat, cmap=cm.gist_heat_r, vmin=vmin, vmax=vmax, extent=[place2[1], place2[2], place1[2], place1[1], ]
                  , rasterized=True)
    ax.set_xticks([]);
    ax.set_yticks([]);
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)

    sns.set_style('ticks')
    newaxs_place2 = add_bigwig_to_axis(place2, vert_bbdict, ax, bins=500, basey=-.35, 
                                       roundby=2, label_for_title="", color=bwcolor,
                                    ylimdict = ylimdict)
    
    newax2 = add_GTF_to_axis(place2, ax, roundby=2, ignore_set=ignore_set, ignore_method=ignore_method, **kwargs)
    # newax2.ticklabel_format(style='sci', axis='x', scilimits=(6, 6))
    s, e = place2[1:]
    for _ in [newax2] + newaxs_place2:
        ticks = [s, (s+e)//2, e]
        _.set_xlim([s, e])
        _.set_xticks(ticks)
    
    # if len(newaxs_place2) > 0:
    #     a = newaxs_place2[0]
    #     print(deepcopy(a.get_xticks()))
    #     newax2.set_xticks(deepcopy(a.get_xticks()))
    
    sns.set_style('ticks')
    newaxs_place1 = add_bigwig_to_axis(place1, vert_bbdict, ax, bins=500, basey=-.68, 
                                       roundby=2, label_for_title="", color=bwcolor,
                                      ylimdict = ylimdict)    
    newax1 = add_GTF_to_L_axis(place1, ax, roundby=2, ignore_set=ignore_set, ignore_method=ignore_method, **kwargs)
    s, e = place1[1:]
    for _ in [newax1] + newaxs_place1:
        ticks = [s, (s+e)//2, e]
        _.set_ylim([s, e])
        _.set_yticks(ticks)
    # newax1.ticklabel_format(style='sci', axis='y', scilimits=(6, 6))
    newax1.set_ylim(*place1[1:])
    newax1.invert_yaxis()
    for a in newaxs_place1:
        a.remove()

    # for a in newaxs_place2:
        # a.ticklabel_format(style='sci', axis='x', scilimits=(6, 6))    

    ax.set_title(name_of_func)
    newax2.set_yticks([])
    newax2.spines['left'].set_visible(False)
    # newax2.set_xticks([])
    return fig, ([newax1, newax2], newaxs_place2, place1, place2), ax
    return ax, ([newax1, newax2], newaxs_place2, place1, place2)


# def make_plots_in_array():
    

def make_scoping_plot_for_fig1_near_diag(coolfile_UPPER, coolfile_LOWER, row, col, get_func_UPPER, get_func_LOWER, name_of_func, vert_bbdict = {},
                                resolution=5000, d = 110000, 
                                d_RR=None, 
                                d_LL=None, 
                                d_RL=None,
                                d_LR=None,
                                useSigma=False, ignore_set=[],
                                vmin=0, vmax=1, balance=True, bwcolor='black', 
                                ylimdict=None, all_ind_to_region=None, 
                                filesoi=None, 
                                ignore_method='ignore',
                                ax = None,
                                **kwargs):
    if d_RR is None:
        d_RR = d
    if d_LL is None:
        d_LL = d
    if d_RL is None:
        d_RL = d
    if d_LR is None:
        d_LR = d

    sns.set(font_scale=1)
    chrom1, s1, e1 = all_ind_to_region[row]
    place1  = (chrom1, s1-d_LL, e1+d_LR)

    chrom2, s2, e2 = all_ind_to_region[col]
    place2  = (chrom2, s2-d_RL, e2+d_RR)

    m_up, _, _ = get_func_UPPER(coolfile_UPPER, place1, place2, filesoi=filesoi, resolution=resolution, balance=balance)
    fullm = np.zeros_like(m_up)
    fullm += np.triu(m_up, k=1)
    
    m_low, _, _ = get_func_LOWER(coolfile_LOWER, place1, place2, filesoi=filesoi, resolution=resolution, balance=balance)
    fullm += np.tril(m_low, k=-1)
    
    

    sns.set_style('ticks')
    if useSigma:
        filtered_mat = scipy.ndimage.gaussian_filter(m, sigma=.75)
    else:
        filtered_mat = fullm
        print(np.max(filtered_mat))
    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 3))
    else:
        fig = None
    if (place2[1] <= place1[2]) and (chrom1 == chrom2) and (place2[2] >= place1[2]):
        ax.matshow(filtered_mat, cmap=cm.gist_heat_r, norm=colors.LogNorm(vmin=vmin, vmax=vmax), extent=[place2[1], place2[2], place1[2], place1[1], ]
                  , rasterized=True)
    else:
        ax.matshow(filtered_mat, cmap=cm.gist_heat_r, vmin=vmin, vmax=vmax, extent=[place2[1], place2[2], place1[2], place1[1], ]
                  , rasterized=True)
    ax.set_xticks([]);
    ax.set_yticks([]);
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)

    sns.set_style('ticks')
    newaxs_place2 = add_bigwig_to_axis(place2, vert_bbdict, ax, bins=500, basey=-.35, 
                                       roundby=2, label_for_title="", color=bwcolor,
                                    ylimdict = ylimdict)
    a = newaxs_place2[0]
    newax2 = add_GTF_to_axis(place2, ax, roundby=2, ignore_set=ignore_set, ignore_method=ignore_method, **kwargs)
    # newax2.ticklabel_format(style='sci', axis='x', scilimits=(6, 6))
    newax2.set_xticks(a.get_xticks())
    newax2.set_xlim(*place2[1:])

    sns.set_style('ticks')
    newaxs_place1 = add_bigwig_to_axis(place1, vert_bbdict, ax, bins=500, basey=-.68, 
                                       roundby=2, label_for_title="", color=bwcolor,
                                      ylimdict = ylimdict)    
    a = newaxs_place1[0]
    newax1 = add_GTF_to_L_axis(place1, ax, roundby=2, ignore_set=ignore_set, ignore_method=ignore_method, **kwargs)
    # newax1.ticklabel_format(style='sci', axis='y', scilimits=(6, 6))
    newax1.set_yticks(deepcopy(a.get_xticks()))
    newax1.set_yticks(deepcopy(a.get_xticks()))
    newax1.set_ylim(*place1[1:])
    newax1.invert_yaxis()
    for a in newaxs_place1:
        a.remove()

    # for a in newaxs_place2:
        # a.ticklabel_format(style='sci', axis='x', scilimits=(6, 6))    

    ax.set_title(name_of_func)

    return fig, ([newax1, newax2], newaxs_place2, place1, place2), ax
    # return ax, ([newax1, newax2], newaxs_place2, place1, place2)


