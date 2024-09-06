import numpy as np
import matplotlib.pyplot as plt
# import pyBigWig
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
import seaborn as sns
from aux_functions import *

mm = 1/2.54/10  # mm in inches

def plot_inds(mat, i1, i2, d, ax=None, trans_func = lambda x: x, **kwargs):
    L1 = max(i1-d, 0)
    L2 = max(i2-d, 0)

    sl1 = slice(L1, i1+d+1)
    sl2 = slice(L2, i2+d+1)

    submat = trans_func(mat[sl1, sl2])
    if ax is None:
        fig, ax = plt.subplots()
    ax.matshow(submat, **kwargs)
    return submat.copy()

def subset_mat(mat, i1, i2, d):
    sl1 = slice(i1-d, i1+d+1)
    sl2 = slice(i2-d, i2+d+1)
    submat = mat[sl1, sl2]
    return submat.copy()

def add_xaxis_labels(label_left, label_right, ax, y=-.1, **kwargs):
    plt.text(0, y, label_left, horizontalalignment='center', verticalalignment='top', transform=ax.transAxes, **kwargs)
    plt.text(1, y, label_right, horizontalalignment='center', verticalalignment='top', transform=ax.transAxes, **kwargs)

def add_yaxis_labels(label_top, label_bot, ax, x = -.1, **kwargs):
    plt.text(x, .95, label_bot, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, **kwargs)
    plt.text(x, .05, label_top, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, **kwargs)


def plot_mat(mat, ax=None, trans_func = lambda x: x, **kwargs):
    submat = trans_func(mat)
    if ax is None:
        fig, ax = plt.subplots()
    ax.matshow(submat, **kwargs)
    return submat


def plot_intrachromosomal_hic(i1, i2, all_ind_to_region, cooldict, ax=None, d=30, plot=True, 
                            balance=True, fgsz=4, useLog=False, n_rows=1, add_chr=False, 
                            ylabel='',
                            xlabel='',
                            res = 50000,
                              **kwargs):
    n = len(cooldict)
    if (ax is None) and (plot==True):
        fig, axs = init_subplots_exact(n, n_rows, fgsz=(fgsz, fgsz), dpi=100)
    l1, l2 = all_ind_to_region[i1], all_ind_to_region[i2]
    if add_chr:
        l1, l2 = add_chr_to_anc(l1), add_chr_to_anc(l2)
    newl1, newl2 = extend_l(l1, d*50000), extend_l(l2, d*50000)
    for ax, (cond, cool) in zip(axs, cooldict.items()):
        try:
            m = cool.matrix(balance=balance).fetch(newl1, newl2)
            if useLog:
                m = np.log(m)
            if plot==True:
                n = len(m)
                _ = plot_mat(m, cmap=cm.gist_heat_r, ax=ax, extent = [n//2, -n//2, -n//2, n//2], **kwargs)
                ax.set_title(f'{cond}: {i1}-{i2}')
        except Exception as e:
            print(e)
            return None
        ax.grid(False)
    for c, a in enumerate(axs):
        n = len(m)
        cutoff = np.round(-n//2*res/1e3/1e3, 1)
        a.set_yticks([-n//2, 0, n//2])
        a.set_xticks([-n//2, 0, n//2])
        a.set_xticklabels([-cutoff, xlabel, cutoff])
        a.tick_params(labeltop = False, top=False)
        a.tick_params(labelbottom = True)
        if c%2 == 0:
            a.set_yticklabels([-cutoff, ylabel, cutoff])
            a.get_yticklabels()[1].set_fontsize(10)
            a.get_yticklabels()[1].set_rotation(90)
            a.set_yticklabels(a.get_yticklabels(), va='center')
        else:
            a.tick_params(labelleft = False)
        a.get_xticklabels()[1].set_fontsize(10)

    return (axs, m, newl1, newl2)


def plot_interchromosomal_pmat_matrices(i1, i2, all_ind_to_region, cooldict_50kb, ax=None, d=30, cool_cond='treg', plot=True, balance=True, **kwargs):
    if (ax is None) and (plot==True):
        fig, ax = plt.subplots()
    l1, l2 = all_ind_to_region[i1], all_ind_to_region[i2]
    newl1, newl2 = extend_l(l1, d*50000), extend_l(l2, d*50000)
    try:
        m = cooldict_50kb[cool_cond].matrix(balance=balance).fetch(newl1, newl2)
        if plot==True:
            _ = plot_mat(m, cmap=cm.gist_heat_r, ax=ax, **kwargs)
            ax.set_title(cool_cond)
            ax.set_title([i1, i2])
    except Exception as e:
        print("Error:", e)
        return None
    return (m, newl1, newl2)

import scipy
import scipy.stats
def make_pval_mat(m1_list, m2_list):
    n = m1_list.shape[2]
    pvals = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            v1 = m1_list[:, i, j]
            v2 = m2_list[:, i, j]
            bad1 = np.isnan(v1); bad2 = np.isnan(v2)
            v1, v2 = v1[~bad1], v2[~bad2]
            pval = scipy.stats.ttest_ind(v1, v2)[1]
            pvals[i, j] = pval
    return pvals

import statsmodels
import statsmodels.nonparametric
import statsmodels.nonparametric.smoothers_lowess
def plot_smoothed_line(y, x, frac=0.66, **kwargs):
    smooth_x, smooth_y = statsmodels.nonparametric.smoothers_lowess.lowess(x, y, frac=frac,).T
    plt.plot(smooth_x, smooth_y, **kwargs)


import matplotlib
import numpy
import pygbrowse

import pygbrowse
import pygbrowse.plots
from pygbrowse import utilities
from pygbrowse.plots import GenomeBrowser
from pygbrowse.plots import *
class MyGenomeBrowser(GenomeBrowser):
    VECTOR_LEGEND_LOC = 0

    def visualize(self, chrom, start, end, ax,
                  fig_width=12,
                  row_heights=1,
                  ax_spacing=0.05,
                  ignore_set = [],
                  minsize=1000, plotfrac=1,
                  num_xticks=10,
                #   seaborn_style=sns.axes_style(style='ticks',
                                                #    rc={'axes.edgecolor': 'w', 'axes.facecolor': '#EAEAF2'}),
                    scilim=True,
                                                   ):
        """
        Generate, display and return a matplotlib.Figure object comprising one or more Axes representing the genomic
        data tracks specified at initialization.
        The region to plot is specified by the parameters chrom, start, and end.
        :param:`fig_width` is specified in inches
        :param:`row_heights` can be specified as a scalar value (in inches), in which case the same row height will be
        used for all subplots, or as an iterable, in which case the row heights will be applied to the subplots
        in order.
        :param chrom:
        :param start:
        :param end:
        :param fig_width:
        :param row_heights:
        :param ax_spacing:
        :param num_xticks:
        :param seaborn_style:
        :return:
        """

        # ToDo: Add gene (or other feature) lookup instead of specifying coordinates.
        start, end = int(start), int(end)

        assert end > start, 'Window end must be greater than window start! Got: {}, {}'.format(start, end)

        # if we receive a scalar here, use it as the height for all rows
        try:
            if len(row_heights) == 1:
                row_heights = row_heights * len(self.subplot_objects)  # treat as a uniform row height
        except TypeError:
            row_heights = [row_heights] * len(self.subplot_objects)  # treat as a uniform row height

        assert len(row_heights) == len(self.subplot_objects)

        span = end - start
        xtick_increment = span / num_xticks
        rounding_increment = 5 * 10 ** np.round(np.log10(xtick_increment) - 1)
        xtick_increment = utilities.roundto(xtick_increment, rounding_increment)
        num_ticks = int(span / xtick_increment) + 1
        round_start = utilities.roundto(start, rounding_increment)

        for ax_idx in range(len(self.subplot_objects)):
            this_ax = ax

            if ax_idx == len(self.subplot_objects) - 1:
                this_ax.set_xticks(np.arange(num_ticks) * xtick_increment + round_start)
                if scilim:
                    this_ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
                # this_ax.set_xlabel('{} position'.format(chrom))
                
            else:  # clear out xticks but plot objects can override this later
                this_ax.set_xlabel('')
                this_ax.set_xticks([])

            plot_object_subset = self.subplot_objects[ax_idx]

            # Set default plot limits (can be changed by client objects)
            this_ax.set_ylim(0, 1)
            this_ax.set_xlim(start, end)
            
            for plot_object in plot_object_subset:
                plot_object.plot(this_ax, chrom=chrom, ws=start, we=end, fig_width=fig_width,
                                 row_height=row_heights[ax_idx], ignore_set=ignore_set, minsize=minsize, plotfrac=plotfrac)
                

            # ToDo: Refactor legend code to get colors and names from objects not from axes handles.
            # if len(this_ax.get_legend_handles_labels()[1]):
            #     this_ax.legend(loc=self.VECTOR_LEGEND_LOC)
        # return fig

class MyGeneModelPlot(pygbrowse.plots.GeneModelPlot):
    def __init__(self, data, ignore_method='grayout', skip_gm=False, skip_rik=False, **kwargs):
        super().__init__(gene_annotation_data = data, **kwargs)
        # self.chevron_width = 0.0004
        # self.feature_height=0.06
        self.chevron_height=0.15/3
        self.chevron_width=0.35/3
        self.chevron_spacing=0.7
        self.ignore_method = ignore_method
        self.skip_gm = skip_gm
        self.skip_rik = skip_rik
        # self.truncation_size=0.10*.001

        # self.utr_endcap_width=0.04*.001*0

    def plot(self, ax, chrom, ws, we, fig_width, row_height, ignore_set=set(), minsize=1000, plotfrac=1, **kwargs):
        # find overlapping genes
        overlapping_genes, overlapping_transcripts, overlapping_components, ids_to_names = self.gene_annotation_data.query(
            chrom, ws, we)
        # overlapping_genes = self.genes.overlapping(chrom, ws, we)
        # overlapping_components = self.components.overlapping(chrom, ws, we)

        gene_display_levels = self._arrange_genes(overlapping_genes.values())
        ax.set_ylim((-0.5, len(gene_display_levels) - 1 + 0.5))

        # convert inches to data coordinates
        chevron_spacing_dt = (we - ws) / (fig_width / self.chevron_spacing)
        chevron_width_dt = (we - ws) / (fig_width / self.chevron_width)
        truncation_width_dt = (we - ws) / (fig_width / self.truncation_size)
        utr_endcap_width_dt = (we - ws) / (fig_width / self.utr_endcap_width)

        feature_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (row_height / self.feature_height)
        chevron_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (row_height / self.chevron_height)
        truncation_height_dt = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (row_height / self.truncation_size)

        for gene_num, level_genes in enumerate(gene_display_levels):

            # ToDo: make this universal. Divide the gene body into non-overlapping segments, each type of which has a template.

            for gene_id in level_genes:
                name_key = 'Name' if 'Name' in overlapping_genes[gene_id] else 'gene_name'
                gene_data = overlapping_genes[gene_id]
                #                 print(gene_id, gene_data['Name'])
                if gene_data[name_key] in ignore_set:
                    ignore=True
                else:
                    ignore=False
                if ('Gm' in gene_data[name_key]) and self.skip_gm:
                    continue
                if ('Rik' in gene_data[name_key]) and self.skip_rik:
                    continue
                    
                # if abs(gene_data['start'] - gene_data['end']) < minsize:
                    # continue
                # if len(self.genes_include) > 0 and (not any(gene_data['Name'] in s for s in self.genes_include)):
                    # continue
                if ignore and (self.ignore_method == 'ignore'):
                    continue
                    # color = 'gray'
                left_truncated = gene_data['start'] < ws
                right_truncated = gene_data['end'] > we

                visible_gene_start = max(gene_data['start'], ws)
                if left_truncated:
                    visible_gene_start += truncation_width_dt * 2
                visible_gene_end = min(gene_data['end'], we)
                if right_truncated:
                    visible_gene_end -= truncation_width_dt * 2

                if ignore and self.ignore_method == 'grayout':
                    plot_color = 'lightgray'
                elif ignore and (self.ignore_method == 'grayout_no_text'):
                    plot_color = 'lightgray'
                else:
                    plot_color = self.color

                ax.plot((visible_gene_start, visible_gene_end), (gene_num, gene_num), c=plot_color, zorder=1)
                if ignore and (self.ignore_method == 'grayout_no_text'):
                    pass
                else:
                    ax.text(x=(visible_gene_start + visible_gene_end) / 2,
                            y=gene_num + feature_height_dt * 1.5,
                            s=gene_data[name_key],
                            ha='center',
                            fontsize=self.gene_name_fontsize)

                num_chevrons = int(max((visible_gene_end - visible_gene_start) / chevron_spacing_dt, 1))
                chevron_remainder = (visible_gene_end - visible_gene_start) - (num_chevrons - 1) * chevron_spacing_dt

                if gene_data['strand'] == '+':
                    chevron_x_delta = -chevron_width_dt
                else:
                    chevron_x_delta = chevron_width_dt

                for chevron_idx in range(num_chevrons):
                    chevron_x = visible_gene_start + chevron_idx * chevron_spacing_dt + chevron_remainder / 2

                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num + chevron_height_dt),
                            color=plot_color)
                    ax.plot((chevron_x, chevron_x + chevron_x_delta), (gene_num, gene_num - chevron_height_dt),
                            color=plot_color)

                # if left_truncated:
                #     y_points = [gene_num, gene_num - truncation_height_dt, gene_num + truncation_height_dt]
                #     left_x_point = ws + 1
                #     right_x_point = ws + truncation_width_dt + 1

                #     x_points = numpy.array([left_x_point, right_x_point, right_x_point])

                #     larr1 = matplotlib.patches.Polygon(numpy.vstack([x_points, y_points]).T,
                #                                        edgecolor=None,
                #                                        facecolor=None,
                #                                        fill=True,
                #                                        transform=ax.transData,
                #                                        zorder=3)
                #     larr2 = matplotlib.patches.Polygon(numpy.vstack([x_points + truncation_width_dt, y_points]).T,
                #                                        edgecolor=None,
                #                                        facecolor=None,
                #                                        fill=True,
                #                                        transform=ax.transData,
                #                                        zorder=3)

                #     ax.add_patch(larr1)
                #     ax.add_patch(larr2)

                # if right_truncated:
                #     y_points = [gene_num, gene_num - truncation_height_dt, gene_num + truncation_height_dt]
                #     left_x_point = we - truncation_width_dt - 1
                #     right_x_point = we - 1

                #     x_points = numpy.array([right_x_point, left_x_point, left_x_point])

                #     rarr1 = matplotlib.patches.Polygon(xy=numpy.vstack([x_points, y_points]).T,
                #                                        edgecolor=None,
                #                                        facecolor='black',
                #                                        fill=True,
                #                                        transform=ax.transData,
                #                                        zorder=3)
                #     rarr2 = matplotlib.patches.Polygon(numpy.vstack([x_points - truncation_width_dt, y_points]).T,
                #                                        edgecolor=None,
                #                                        facecolor='black',
                #                                        fill=True,
                #                                        transform=ax.transData,
                #                                        zorder=3)
                #     ax.add_patch(rarr1)
                #     ax.add_patch(rarr2)

                # # Identify components belonging to this gene
                # # this_gene_components = set([])
                # # for transcript_id in gene_data['transcripts']:
                # #     for component_id in overlapping_components:
                # #         this_gene_components.add(component_id)

                # plot components
                for component_id in overlapping_components:
                    component_data = overlapping_components[component_id]
                
                    if ((component_data['start'] >= visible_gene_start) and (
                            component_data['start'] <= visible_gene_end)) or (
                            (component_data['end'] >= visible_gene_start) and (
                            component_data['end'] <= visible_gene_end)):

                        # ToDo: systematize and condense the following:
                        if component_data['type'] == 'five_prime_UTR':
                            # plot the "body" of the UTR
                            if gene_data['strand'] == '+':
                                utr_body = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt / 2),
                                    width=max(component_data['end'] - component_data['start'], 0)*10,
                                    height=feature_height_dt,
                                    facecolor=plot_color,
                                    linewidth=0)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['end'] - utr_endcap_width_dt, gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=plot_color,
                                    linewidth=0)

                            else:

                                utr_body = matplotlib.patches.Rectangle(xy=(
                                    component_data['start'] + utr_endcap_width_dt, gene_num - feature_height_dt / 2),
                                    width=max(
                                        component_data['end'] - component_data[
                                            'start'] - utr_endcap_width_dt, 0),
                                    height=feature_height_dt,
                                    facecolor=self.color,
                                    linewidth=0)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=plot_color,
                                    linewidth=0)

                            ax.add_patch(utr_body)
                            ax.add_patch(utr_endcap)

                        elif component_data['type'] == 'three_prime_UTR':
                            # plot the "body" of the UTR
                            if gene_data['strand'] == '-':
                                utr_body = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt / 2),
                                    width=max(component_data['end'] - component_data['start'] - utr_endcap_width_dt, 0),
                                    height=feature_height_dt,
                                    facecolor=plot_color, linewidth=0)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['end'] - utr_endcap_width_dt, gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=plot_color, linewidth=0)

                            else:
                                utr_body = matplotlib.patches.Rectangle(xy=(
                                    component_data['start'] + utr_endcap_width_dt, gene_num - feature_height_dt / 2),
                                    width=max(
                                        component_data['end'] - component_data[
                                            'start'] - self.utr_endcap_width, 0),
                                    height=feature_height_dt,
                                    facecolor=plot_color, linewidth=0)
                                utr_endcap = matplotlib.patches.Rectangle(
                                    xy=(component_data['start'], gene_num - feature_height_dt),
                                    width=utr_endcap_width_dt,
                                    height=feature_height_dt * 2,
                                    facecolor=plot_color, linewidth=0)

                            ax.add_patch(utr_body)
                            ax.add_patch(utr_endcap)

                        elif component_data['type'] == 'CDS':
                            cds = matplotlib.patches.Rectangle(
                                xy=(component_data['start'], gene_num - feature_height_dt),
                                width=component_data['end'] - component_data['start'],
                                height=feature_height_dt * 2,
                                facecolor=plot_color, linewidth=0)
                            ax.add_patch(cds)

        ax.set_yticks([])
        ax.set_ylabel(self.label)

def add_GTF_to_axis(place, ax, ignore_set=[], roundby=-5, gene_name_fontsize=14, outpath= '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz', 
                    bound_L = -.15, bound_R = .1, ignore_method='ignore', skip_rik=False, skip_gm=False, **kwargs):
    chrom1, s1, e1 = place
    genemodels = pygbrowse.datasources.Gff3Annotations(f'{outpath}.bgzf')
    genemodel_plotter = MyGeneModelPlot(genemodels, gene_name_fontsize=gene_name_fontsize, ignore_method=ignore_method,
                                        skip_rik=skip_rik, skip_gm=skip_gm)
    genes_only_plotter = MyGenomeBrowser([[genemodel_plotter]])

    tmpax = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, -.15, 1, .1))
    _ = genes_only_plotter.visualize(f'chr{chrom1}', s1, e1, ax=tmpax, ignore_set = ignore_set, **kwargs)
    
    s, e = tmpax.get_xlim()
    midpt = (e+s)//2
    tmpax.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])
    tmpax.set_ylim(-.25, .75)
    tmpax.spines['top'].set_visible(False)
    tmpax.spines['right'].set_visible(False)

    newax1 = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, bound_L, 1, bound_R))

    for a in [tmpax]:
        s, e = a.get_xlim()
        midpt = (e+s)//2
        a.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])
    for j in tmpax.lines:
        x, y = j.get_xdata(), j.get_ydata()
        newax1.plot(x, y, label=j.get_label(), linestyle=j.get_linestyle(), color=j.get_color(), linewidth=j.get_linewidth(),
                    zorder=1)

    for patch in tmpax.patches:
        if type(patch) == mpl.patches.Rectangle:
            x, y = patch.get_xy()
            height = patch.get_height()
            width = patch.get_width()
            patch = mpl.patches.Rectangle((x, y), width, height, edgecolor=patch.get_edgecolor(), facecolor=patch.get_facecolor(),
                                          fill=patch.get_fill(), linewidth=patch.get_linewidth(), linestyle=patch.get_linestyle(),
                                          alpha=patch.get_alpha(), zorder=2)
            newax1.add_patch(patch)
        else:
            print("Not writing", type(patch))
            continue    
    tmpax.remove()
    newax1.ticklabel_format(axis='both', style='sci', scilimits=(6, 6))
    newax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    newax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    newax1.set_yticks([])
    newax1.grid(False)

    newax1.spines['top'].set_visible(False)
    newax1.spines['right'].set_visible(False)
    for child in tmpax.get_children():
        if type(child) == mpl.text.Text:
            x, y, text = child._x, child._y, child.get_text()
            if len(text) > 1:
                newax1.text(x, y, text, zorder=2, fontsize=gene_name_fontsize, rotation=0)
    newax1.set_ylim(-.5, 1)
    newax1.set_xlabel(f"Mb, chr{chrom1}")

    return newax1

import matplotlib as mpl
from copy import deepcopy
import seaborn as sns
def add_GTF_to_L_axis(place, ax, ignore_set=[], roundby=-5, gene_name_fontsize=14, 
                      bound_L = -.15, bound_R = .1, skip_rik=False, skip_gm=False,
                      outpath= '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz', ignore_method = 'ignore', **kwargs):
    chrom1, s1, e1 = place
    genemodels = pygbrowse.datasources.Gff3Annotations(f'{outpath}.bgzf')
    genemodel_plotter = MyGeneModelPlot(genemodels, gene_name_fontsize=gene_name_fontsize, ignore_method=ignore_method,
                                        skip_rik=skip_rik, skip_gm=skip_gm,)
    genes_only_plotter = MyGenomeBrowser([[genemodel_plotter]], )

    newax1 = ax.inset_axes(transform=ax.transAxes,
        bounds = (bound_L, 0, bound_R, 1))
    newax1.invert_yaxis()

    tmpax = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, bound_L, 1, bound_R))

    _ = genes_only_plotter.visualize(f'chr{chrom1}', s1, e1, ax=tmpax, ignore_set = ignore_set, **kwargs)
    
    for a in [tmpax]:
        s, e = a.get_xlim()
        midpt = (e+s)//2
        a.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])

    for j in tmpax.lines:
        x, y = j.get_ydata(), j.get_xdata()
        newax1.plot(x, y, label=j.get_label(), linestyle=j.get_linestyle(), color=j.get_color(), linewidth=j.get_linewidth(),
                    zorder=1)

    for patch in tmpax.patches:
        if type(patch) == mpl.patches.Rectangle:
            y, x = patch.get_xy()
            width = patch.get_height()
            height = patch.get_width()
            patch = mpl.patches.Rectangle((x, y), width, height, edgecolor=patch.get_edgecolor(), facecolor=patch.get_facecolor(),
                                          fill=patch.get_fill(), linewidth=patch.get_linewidth(), linestyle=patch.get_linestyle(),
                                          alpha=patch.get_alpha(), zorder=2)
            newax1.add_patch(patch)
        else:
            print("Not writing", type(patch))
            continue    


    newax1.spines['top'].set_visible(False)
    newax1.spines['right'].set_visible(False)
            
    for child in tmpax.get_children():
        if type(child) == mpl.text.Text:
            x, y, text = child._x, child._y, child.get_text()
            if len(text) > 1:
                newax1.text(y, x, text, zorder=2, fontsize=gene_name_fontsize, rotation=-90)
    
    newax1.set_xlim(tmpax.get_ylim())
    newax1.set_ylim(tmpax.get_xlim())

    newax1.set_xlabel(deepcopy(tmpax.get_ylabel()))
    newax1.set_ylabel(deepcopy(tmpax.get_xlabel()))
    newax1.set_xticks(deepcopy(tmpax.get_yticks()))
    newax1.set_yticks(deepcopy(tmpax.get_xticks()))

    newax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    newax1.invert_yaxis()

    # newax1.set_xlim([-1, 1])
    tmpax.remove()
    return newax1


def add_one_GTF_to_axis(place1, ax, ignore_set=[], roundby=-5, gene_name_fontsize=12):
    chrom1, s1, e1 = place1
    outpath = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz'
    genemodels = pygbrowse.datasources.Gff3Annotations(f'{outpath}.bgzf')
    genemodel_plotter = pygbrowse.plots.GeneModelPlot(genemodels, gene_name_fontsize=gene_name_fontsize)
    genes_only_plotter = MyGenomeBrowser([[genemodel_plotter]])

    newax1 = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, -.3, 1, .2))
    _ = genes_only_plotter.visualize(f'chr{chrom1}', s1, e1, ax=newax1, ignore_set = ignore_set)
    for a in [newax1]:
        s, e = a.get_xlim()
        midpt = (e+s)//2
        a.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])
    return newax1


def plot_one_GTF_on_axis(place, ax, ignore_set=[], roundby=-5, gene_name_fontsize=14, outpath= '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz', 
                    ignore_method='ignore', skip_gm = False, skip_rik = False, **kwargs):
    chrom1, s1, e1 = place
    genemodels = pygbrowse.datasources.Gff3Annotations(f'{outpath}.bgzf')
    genemodel_plotter = MyGeneModelPlot(genemodels, gene_name_fontsize=gene_name_fontsize, ignore_method=ignore_method,
                                        skip_gm = skip_gm, skip_rik = skip_rik)
    genes_only_plotter = MyGenomeBrowser([[genemodel_plotter]])

    tmpax = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, 0, 1, 1))
    _ = genes_only_plotter.visualize(f'chr{chrom1}', s1, e1, ax=tmpax, ignore_set = ignore_set, **kwargs)
    
    s, e = tmpax.get_xlim()
    midpt = (e+s)//2
    tmpax.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])
    tmpax.set_ylim(-.25, .75)
    tmpax.spines['top'].set_visible(False)
    tmpax.spines['right'].set_visible(False)

    newax1 = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, 0, 1, 1))

    for a in [tmpax]:
        s, e = a.get_xlim()
        midpt = (e+s)//2
        a.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])
    for j in tmpax.lines:
        x, y = j.get_xdata(), j.get_ydata()
        newax1.plot(x, y, label=j.get_label(), linestyle=j.get_linestyle(), color=j.get_color(), linewidth=j.get_linewidth(),
                    zorder=1)

    for patch in tmpax.patches:
        if type(patch) == mpl.patches.Rectangle:
            x, y = patch.get_xy()
            height = patch.get_height()
            width = patch.get_width()
            patch = mpl.patches.Rectangle((x, y), width, height, edgecolor=patch.get_edgecolor(), facecolor=patch.get_facecolor(),
                                          fill=patch.get_fill(), linewidth=patch.get_linewidth(), linestyle=patch.get_linestyle(),
                                          alpha=patch.get_alpha(), zorder=2)
            newax1.add_patch(patch)
        else:
            print("Not writing", type(patch))
            continue    
    tmpax.remove()
    newax1.ticklabel_format(axis='both', style='sci', scilimits=(6, 6))
    newax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    newax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    newax1.set_yticks([])
    newax1.grid(False)

    newax1.spines['top'].set_visible(False)
    newax1.spines['right'].set_visible(False)
    for child in tmpax.get_children():
        if type(child) == mpl.text.Text:
            x, y, text = child._x, child._y, child.get_text()
            if len(text) > 1:
                newax1.text(x, y, text, zorder=2, fontsize=gene_name_fontsize, rotation=0)
    newax1.set_ylim(-.5, 1)
    newax1.set_xlabel(f"Mb, chr{chrom1}")

    return newax1



# def plot_one_GTF_on_axis(place1, ax, ignore_set=[], roundby=-5, gene_name_fontsize=12, outpath = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz', 
#                          skip_gm = False, skip_rik = False, **kwargs):
#     chrom1, s1, e1 = place1
#     genemodels = pygbrowse.datasources.Gff3Annotations(f'{outpath}.bgzf')
#     genemodel_plotter = MyGeneModelPlot(genemodels, gene_name_fontsize=gene_name_fontsize,
#                                         skip_gm = skip_gm, skip_rik = skip_rik)
#     genes_only_plotter = MyGenomeBrowser([[genemodel_plotter]])

#     tmpax = ax.inset_axes(transform=ax.transAxes,
#         bounds = (0, 0, 1, 1)); 
#     _ = genes_only_plotter.visualize(f'{chrom1}', s1, e1, ax=tmpax, ignore_set = ignore_set, **kwargs)
    
#     s, e = tmpax.get_xlim()
#     midpt = (e+s)//2
#     tmpax.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])
#     tmpax.set_ylim(-.25, .75)
#     tmpax.spines['top'].set_visible(False)
#     tmpax.spines['right'].set_visible(False)

#     newax1 = ax

#     for a in [tmpax]:
#         s, e = a.get_xlim()
#         midpt = (e+s)//2
#         a.set_xticks([np.round(s, roundby), np.round(midpt, roundby), np.round(e, roundby)])
#     for j in tmpax.lines:
#         newax1.plot(j.get_xdata(), j.get_ydata(), color='black', zorder=-1)
#     for patch in tmpax.patches:
#         if type(patch) == mpl.patches.Rectangle:
#             x, y = patch.get_xy()
#             height = patch.get_height()
#             width = patch.get_width()
#             patch = mpl.patches.Rectangle((x, y), width, height, color='black', zorder=2)
#             newax1.add_patch(patch)
#         else:
#             print("Not writing", type(patch))
#             continue    
#     tmpax.remove()
#     newax1.spines['top'].set_visible(False)
#     newax1.spines['right'].set_visible(False)
#     for child in tmpax.get_children():
#         if type(child) == mpl.text.Text:
#             x, y, text = child._x, child._y, child.get_text()
#             if len(text) > 1:
#                 newax1.text(x, y, text, zorder=2, fontsize=gene_name_fontsize, rotation=0)
#                 print(text)
#     newax1.set_ylim(-.5, 1)
#     return newax1


from aux_functions import make_int
def add_motifs_to_ax(newax, place1, motif_bedtool):
    chrom1, s1, e1 = make_int(place1)
    for row in motif_bedtool:
        chrom, s, e = make_int(row)[:3]
        pval = -np.log10(float(row[4]))
        alpha = (pval-4)/6
        newax.axvspan(s, e, alpha=alpha, color='blue')
    return newax



def add_chr_to_anc(anc):
    anc = list(anc)
    if 'chr' not in anc[0]:
        anc[0] = 'chr' + anc[0]
    return tuple(anc)

def add_bigwig_to_axis(place, bbdict, ax, ignore_set=[], bins=1000, basey=-1, 
                     roundby = -5, label_for_title='', ylimdict = {}, color=None, **kwargs):
    newaxs = []
    for c, (label, bb) in enumerate(bbdict.items()):
        y = basey - .15*c
        newax1 = ax.inset_axes(transform=ax.transAxes,
            bounds = (0, y, 1, .1))
        xs = np.linspace(place[1], place[2], bins)
        if 'chr' in list(bb.chromsizes)[0]:
            v = bb.fetch(*add_chr_to_anc(place), bins=bins)
        else:
            v = bb.fetch(*place, bins=bins)
        newax1.plot(xs, v, color=color, **kwargs)
        newax1.set_xlim(xs[0], xs[-1])
        
        newax1.set_ylim(ylimdict.get(label, [0, 80]))
        s, e = newax1.get_ylim()
        midpt = (e+s)/2
        newax1.set_yticks([s, midpt, e])
        
        s, e = newax1.get_xlim()
        midpt = (e+s)/2
        newax1.set_xlim(s, e)
        newax1.set_ylabel(label_for_title + label, rotation=0, ha='left', y=.65)
        newax1.yaxis.set_label_position("right")
        newax1.spines['top'].set_visible(False)
        newax1.spines['right'].set_visible(False)
        newax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
        newax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
        newaxs.append(newax1)

    return newaxs

def add_vector_to_L_axis(vector, ax, ignore_set=[], bins=1000, x=-1, 
                     roundby = -5, label_for_title='', ylimdict = {}, cmap=None, **kwargs):
    newaxs = []
    newax1 = ax.inset_axes(transform=ax.transAxes,
        bounds = (x, 0, .1, 1))
    newax1.matshow(vector[:, None], cmap=cmap, **kwargs)


def add_vector_to_axis(vector, ax, ignore_set=[], bins=1000, basey=-1, 
                     roundby = -5, label_for_title='', ylimdict = {}, cmap=None, **kwargs):
    newaxs = []

    y = basey
    newax1 = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, y, 1, .1))
    newax1.matshow(vector[:, None], cmap=cmap, **kwargs)
    # newax1.set_xlim(xs[0], xs[-1])
    
    # newax1.set_ylim(ylimdict.get(label, [0, 80]))
    # s, e = newax1.get_ylim()
    # midpt = (e+s)/2
    # newax1.set_yticks([s, midpt, e])
    
    # s, e = newax1.get_xlim()
    # midpt = (e+s)/2
    # newax1.set_xlim(s, e)
    # newax1.set_ylabel(label_for_title + label, rotation=0, ha='left', y=.65)
    # newax1.yaxis.set_label_position("right")
    # newax1.spines['top'].set_visible(False)
    # newax1.spines['right'].set_visible(False)
    # newax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    # newax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    # newaxs.append(newax1)

    # return newaxs



def add_bigwig_to_axis_fill(place, bbdict, ax, ignore_set=[], bins=1000, basey=-1, 
                     roundby = -5, label_for_title='', ylimdict = {}, color=None, use_log=False, delta=.15, **kwargs):
    newaxs = []
    for c, (label, bb) in enumerate(bbdict.items()):
        y = basey - delta*c
        newax1 = ax.inset_axes(transform=ax.transAxes, sharex=ax,
            bounds = (0, y, 1, (delta*2)/3))
        xs = np.linspace(place[1], place[2], bins)
        if 'chr' in list(bb.chromsizes)[0]:
            v = bb.fetch(*add_chr_to_anc(place), bins=bins)
        else:
            v = bb.fetch(*place, bins=bins)
        if use_log:
            v = np.log10(v+1)
        if isinstance(color, list):
            newax1.fill_between(xs, v, color=color[c], **kwargs)
        else:
            newax1.fill_between(xs, v, color=color, **kwargs)
        newax1.set_xlim(xs[0], xs[-1])
        
        newax1.set_ylim(ylimdict.get(label, [0, 80]))
        s, e = newax1.get_ylim()
        midpt = (e+s)/2
        newax1.set_yticks([s, midpt, e])
        
        s, e = newax1.get_xlim()
        midpt = (e+s)/2
        newax1.set_xlim(s, e)
        newax1.set_ylabel(label_for_title + label, rotation=0, ha='left', y=.65)
        newax1.yaxis.set_label_position("right")
        newax1.spines['top'].set_visible(False)
        newax1.spines['right'].set_visible(False)
        newaxs.append(newax1)
    return newaxs

def plot_bigwig_on_axis(place, bbdict, newax1, condition, ignore_set=[], 
                        bins=1000, basey=-1, roundby = -5, label_for_title='', ylimdict = {},
                        color='blue', alpha=.2, default_ymax = 150, sizefactor=1):
    if isinstance(bbdict, dict):
        bb = bbdict[condition]
    else:
        bb = bbdict

    xs = np.linspace(place[1], place[2], bins)
    if 'chr' in list(bb.chromsizes)[0]:
        v = bb.fetch(*add_chr_to_anc(place), bins=bins)
    else:
        v = bb.fetch(*place, bins=bins)
    v = np.ravel(v)/sizefactor
    newax1.plot(xs, v, color=color, alpha=alpha, linewidth=1)
    newax1.set_xlim(xs[0], xs[-1])
    newax1.set_ylim(ylimdict.get(condition, [0, default_ymax]))
    s, e = newax1.get_ylim()
    midpt = (e+s)/2
    newax1.set_yticks([e])
    s, e = newax1.get_xlim()
    midpt = (e+s)/2
    xticks = np.linspace(s, e, 10)
    newax1.set_xticks(xticks)
    newax1.set_xlim(s, e)
    newax1.set_ylabel(label_for_title, rotation=0, ha='left')
    newax1.tick_params(axis='y', labelright=False, right=False, labelleft=True, left=True)
    newax1.yaxis.set_label_position("right")

def plot_bigwig_on_axis_fill(place, bbdict, newax1, condition, ignore_set=[], 
                        bins=1000, basey=-1, roundby = -5, label_for_title='', ylimdict = {},
                        color='blue', alpha=.2, default_ymax = 150, chrom_converter={}):
    if isinstance(bbdict, dict):
        bb = bbdict[condition]
    else:
        bb = bbdict
    xs = np.linspace(place[1], place[2], bins)
    if any(['chr' in x for x in  list(bb.chromsizes)]):
        v = bb.fetch(*add_chr_to_anc(place), bins=bins)
    elif place[0] not in list(bb.chromsizes.keys()):
        newplace = list(place)
        newplace[0] = chrom_converter.get(place[0], place[0])
        newplace = tuple(newplace)
        v = bb.fetch(*newplace, bins=bins)
    else:
        v = bb.fetch(*place, bins=bins)
    newax1.fill_between(xs, v, color=color)
    newax1.set_xlim(xs[0], xs[-1])
    newax1.set_ylim(ylimdict.get(condition, [0, default_ymax]))
    s, e = newax1.get_ylim()
    midpt = (e+s)/2
    newax1.set_yticks([e])
    s, e = newax1.get_xlim()
    midpt = (e+s)/2
    xticks = np.linspace(s, e, 10)
    newax1.set_xticks(xticks)
    newax1.set_xlim(s, e)
    newax1.set_ylabel(label_for_title, rotation=0, ha='left')
    newax1.tick_params(axis='y', labelright=False, right=False, labelleft=True, left=True)
    newax1.yaxis.set_label_position("right")


def add_arrows_to_ax(rows, cols,  ax):
    shift = 15
    for c, i in enumerate(zip(rows, cols)):
        tmprow, tmpcol = i
        if c == 0: label = 'Megaloop foci'; 
        else: label='';
        ax.arrow(tmpcol-shift, tmprow, shift*.65, 0, color='black', 
                head_width=shift//4, length_includes_head=True,
                width=1, label=label, )
    ax.legend()



def add_loop_anchors_to_ax(place1, place2, anchors, ax, res=5000):
    chrom1, s1, e1 = place1
    chrom2, s2, e2 = place2
    delta = e2-s1
    loopax_x = ax.inset_axes(transform=ax.transAxes,
        bounds = (0, -.06, 1, .05)); loopax_x.set_xlim(ax.get_xlim())

    anc_xs = []
    for i in anchors.intersect(pbt.BedTool([[chrom1, s1, e1]]), u=True):
        anc_xs.append((int(i[1])-s1)/res)
    loopax_x.vlines(anc_xs, -10, 0, color='black')
    loopax_x.set_ylabel("Short-scale loops", rotation=0, va='bottom', ha='left')
    loopax_x.yaxis.set_label_position(position='right')

    loopax_y = ax.inset_axes(transform=ax.transAxes,
        bounds = (-.06, 0, .05, 1)); loopax_y.set_ylim(ax.get_ylim())

    anc_xs = []
    for i in anchors.intersect(pbt.BedTool([[chrom2, s2, e2]]), u=True):
        anc_xs.append((int(i[1])-s2)/res)
    loopax_y.hlines(anc_xs, -10, 0, color='black')

    for a in [loopax_x, loopax_y]:
        a.set_yticks([])
        a.set_xticks([])
        for spine in a.spines:
            a.spines[spine].set_visible(True)


def add_loop_anchors_to_ax_AXVSPAN(place1, anchors, ax):
    chrom1, s1, e1 = place1
    for c, i in enumerate(anchors.intersect(pbt.BedTool([place1]), u=True)):
        s, e = i[1], i[2]
        s, e = map(int, [s, e])
        if c == 0:
            ax.axvspan(s, e, color='black', alpha=.15, linewidth=0, label='Short-scale loop anchors')
        else:
            ax.axvspan(s, e, color='black', alpha=.15, linewidth=0)
    return ax

def add_loop_anchors_to_ax_AXVSPAN_LEFT(place1, anchors, ax, res=5000):
    chrom1, s1, e1 = place1
    for i in anchors.intersect(pbt.BedTool([place1]), u=True):
        s, e = i[1], i[2]
        s, e = map(int, [s, e])
        ax.axhspan(s, e, color='black', alpha=.15, linewidth=0, label='Short-scale loop anchors')
    return ax


def add_pearson_to_plot(X, Y, ax, textcolor = 'black', yloc = .98, label='', test=scipy.stats.ks_2samp, **kwargs):
    # bad = np.isnan(X) | np.isnan(Y) | np.isinf(X) | np.isinf(Y)
    r, pval = test(X, Y)
    ax.text(0.02, yloc, f'r={r.round(2)}\np={format_pvalue(pval)}', transform = ax.transAxes, va='top', fontsize=6, color=textcolor, ha='left')

def add_pval_to_plot(X, Y, ax, textcolor = 'black', xloc = .02, yloc = .98, label='', test=scipy.stats.ks_2samp, fontsize=6, **kwargs):
    # bad = np.isnan(X) | np.isnan(Y) | np.isinf(X) | np.isinf(Y)
    r, pval = test(X, Y)
    ax.text(xloc, yloc, f'{label} p={format_pvalue(pval)}', transform = ax.transAxes, va='top', fontsize=fontsize, color=textcolor)

def scatter_with_pearson(X, Y, ax, textcolor = 'black', yloc = .98, xloc=.02, fontsize=6, test=scipy.stats.pearsonr,
                         **kwargs):
    bad = np.isnan(X) | np.isnan(Y) | np.isinf(X) | np.isinf(Y)
    r, pval = test(X[~bad], Y[~bad])

    ax.text(xloc, yloc, f'r={r.round(2)}\np={format_pvalue(pval)}', transform = ax.transAxes, va='top', fontsize=fontsize, color=textcolor)
    c = ax.scatter(X, Y, **kwargs)
    return r, pval, c

def plot_baseline_chip_values(geneset_dict, color_dict, values_at_all_sites, metadata_at_all_sites, conds, peaktype = 'ATAC', min_number = 10, legend=True):
    figs = []
    for title, cond1 in conds.items():
        plt.figure(figsize=(6, 4))
        values1 = values_at_all_sites[cond1]
        value_meta1 = metadata_at_all_sites[cond1]
        xs_treg = values1
        baseline_vec = scipy.ndimage.gaussian_filter1d(np.log2(xs_treg+1), sigma=10)
        baseline = baseline_vec.mean(axis=0)
        xs = np.linspace(-2500, 2500, len(baseline))
        plt.plot(xs, baseline, label='All peaks', color = 'black', zorder=1)

        for name, genes in geneset_dict.items():
            v = values1[value_meta1['gene'].isin(genes)]
            xs_treg = scipy.ndimage.gaussian_filter1d(np.log2(v+1), sigma=10)
            n_x = len(value_meta1['gene'][value_meta1['gene'].isin(genes)].unique())
            if n_x < min_number:
                continue
            x = xs_treg.mean(axis=0)
            xs = np.linspace(-2500, 2500, len(x))
            pval = scipy.stats.ranksums(xs_treg[:, 500], baseline_vec[:, 500])[1]
            plt.plot(xs, x, label=f'{name}; \nn={n_x}; p={pval:.2e}', color = color_dict[name], zorder=2)
            plt.xlabel(f"{peaktype} peaks near DEGS")
            plt.ylabel("Signal")
        figs.append(plt.gcf())
        plt.title(title)    
        if legend:
            plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.tight_layout()
    return figs

def format_pvalue(pval):
    if pval > 0.05:
        return 'NS'
    elif pval > .001:
        return f'{pval:.3f}'
    else:
        return f'{pval:.1e}'

def plot_comparison_chip_values(geneset_dict, color_dict, values_at_all_sites, metadata_at_all_sites, conds, peaktype='ATAC', use_log=True, ax=None):
    vals = {} 
    for title, (cond1, cond2) in conds.items():
        if ax is None:
            plt.figure(figsize=(8, 4))
        else:
            plt.sca(ax)
        values1 = values_at_all_sites[cond1]
        values2 = values_at_all_sites[cond2]
        value_meta1 = metadata_at_all_sites[cond1]
        value_meta2 = metadata_at_all_sites[cond2]

        if use_log:
            xs_treg = np.log2(values1+1)
            xs_tcon = np.log2(values2+1)
        else:
            xs_treg = values1
            xs_tcon = values2
        baseline_vec = xs_treg-xs_tcon
        baseline = scipy.ndimage.gaussian_filter1d(xs_treg.mean(axis=0) - xs_tcon.mean(axis=0), sigma=10)
        xs = np.linspace(-2500, 2500, len(baseline))
        plt.plot(xs, baseline, label='All peaks', color = 'black', zorder=1)
        for name, genes in geneset_dict.items():
            xs_treg = np.log2(values1[value_meta1['gene'].isin(genes)]+1)
            xs_tcon = np.log2(values2[value_meta2['gene'].isin(genes)]+1)
            n_x = len(value_meta1['gene'][value_meta1['gene'].isin(genes)].unique())
            if n_x < 10:
                continue
            x = scipy.ndimage.gaussian_filter1d(xs_treg.mean(axis=0) - xs_tcon.mean(axis=0), sigma=10)
            xs = np.linspace(-2500, 2500, len(x))
            n = len(genes)
            pval = scipy.stats.ranksums((xs_treg-xs_tcon)[:, 500], baseline_vec[:, 500])[1]
            pval = format_pvalue(pval)
            plt.plot(xs, x, label=f'{name}; n={n}; \n p={pval}', color = color_dict[name], zorder=2, )

        plt.xlabel(f"{peaktype} peaks near DEGS")
        plt.ylabel("Signal")
        plt.title(f"{title} signal \nat {peaktype} peaks near DEGs")
        plt.legend(bbox_to_anchor=(1, 1), loc = 'upper left', fontsize=10) 
        plt.tight_layout()
    
def loop_grange_to_granges(loop_grange):
    return loop_grange.split('|')

def loop_grange_to_tuple(loop_grange):
    l1, l2 = loop_grange_to_granges(loop_grange)
    l1, l2 = map(grange_to_tuple, [l1, l2])
    l1, l2 = map(list, [l1, l2])
    l = l1 + l2
    return l


from aux_functions import extend_l
import cooler
from aux_functions import remove_chr, add_chr_to_anc

def fetch_from_cooler(cool, l1, l2, balance=True):
    if 'chr' in cool.chromnames[0]:
        l1, l2 = map(add_chr_to_anc, [l1, l2])
    else:
        l1, l2 = map(remove_chr, [l1, l2])
    m = cool.matrix(balance=balance).fetch(l1, l2)
    return m

def intersects_diagonal(l1, l2):
    if (l2[1] < l1[2]) and (l2[2] > l1[1]):
        return True
    else:
        return False

import math
from scipy.ndimage import gaussian_filter

def plot_from_multiple_coolers_and_one_grange(cooldict, grange, d, axs = None, loopdict=None, 
                                                useSigma=True, sigma=1, balance=True, ylabel="", xlabel=""):
    if type(cooldict) != dict:
        cooldict = {'cool' : cooldict}
    n = len(cooldict)
    rows = 2
    cols = math.ceil(n/rows)
    if axs is None:
        fig, axs = init_subplots_exact(n, 1, fgsz=(40*mm, 40*mm), dpi = 150)
    axs = np.ravel(axs)
    if n == 1:
        axs = [axs]
    for c, (name, coolfile) in enumerate(cooldict.items()):
        l1, l2 = loop_grange_to_granges(grange)
        l1, l2 = map(grange_to_tuple, [l1, l2])
        l1, l2 = map(make_int, [l1, l2])
        l1, l2 =  extend_l(l1, d), extend_l(l2, d)
        m = fetch_from_cooler(coolfile, l1, l2, balance=balance).astype(float)
        if useSigma:
            without_nan = m
            without_nan[np.isnan(m)] = 0.0
            nan_indicator = (~np.isnan(m)).astype(float)
            without_nan = gaussian_filter(without_nan, sigma=sigma)
            nan_indicator = gaussian_filter(nan_indicator, sigma=sigma)
            m = without_nan/nan_indicator
        if intersects_diagonal(l1, l2):
            axs[c].matshow(m, cmap='gist_heat_r', norm=colors.LogNorm(), extent=[l2[1], l2[2], l1[2], l1[1]])
        else:
            axs[c].matshow(m, cmap='gist_heat_r', extent=[l2[1], l2[2], l1[2], l1[1]])

        axs[c].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
        axs[c].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
        axs[c].set_xticks([l2[1], l2[2]])
        axs[c].set_yticks([l1[1], l1[2]])
        plt.sca(axs[c])
        plt.text(1, -.1, 'chr' + l2[0], ha='right', va='top', transform=axs[c].transAxes, fontsize=6)
        if c == 0:
            axs[c].set_ylabel(ylabel)
            plt.text(-.2, 0, 'chr' + l1[0], ha='right', va='top', transform=axs[c].transAxes, fontsize=6)
        else:
            axs[c].tick_params(labelleft=False, left=False)

        axs[c].tick_params(labeltop=False, top=False, bottom=True, labelbottom=True)
        axs[c].set_xlabel(xlabel)
        axs[c].grid(False)
        plt.title(name)
        if loopdict is not None:
            decorate_ax_with_loops(axs[c], l1, l2, loopdict[name])
    return fig

import math
def plot_from_one_cooler_and_multiple_granges(cooldict, granges, d, axs = None, loopdict=None, cond=None, cols=6, suptitle=None,
                                              add_chr=False):
    if type(cooldict) != dict:
        cooldict = {'cool' : cooldict}
    coolfile = cooldict[cond]
    n = len(granges)
    if axs is None:
        fig, axs = init_subplots_exact(n, 1, fgsz=(40*mm, 40*mm), dpi = 200)
        axs = np.ravel(axs)
    for c, (grange) in enumerate(granges):
        l1, l2 = loop_grange_to_granges(grange)
        l1, l2 = map(grange_to_tuple, [l1, l2])
        l1, l2 = map(make_int, [l1, l2])
        l1, l2 =  extend_l(l1, d), extend_l(l2, d)
        m = fetch_from_cooler(coolfile, l1, l2)
        if intersects_diagonal(l1, l2):
            axs[c].matshow(m, cmap='gist_heat_r', norm=colors.LogNorm(), extent=[l2[1], l2[2], l1[2], l1[1]])
        else:
            axs[c].matshow(m, cmap='gist_heat_r', extent=[l2[1], l2[2], l1[2], l1[1]])
        if loopdict is not None:
            decorate_ax_with_loops(axs[c], l1, l2, loopdict[cond], s = 50, linewidth=.5)
        if add_chr:
            axs[c].set_xlabel(f"chr{l1[0]}:{l1[1]//1000:,} - {l1[2]//1000:,} (kb)")
            axs[c].set_ylabel(f"chr{l2[0]}:{l2[1]//1000:,} - {l2[2]//1000:,} (kb)")
        else:
            axs[c].set_xlabel(f"{l1[0]} : {l1[1]} - {l1[2]}")
            axs[c].set_ylabel(f"{l2[0]} : {l2[1]} - {l2[2]}")
        axs[c].set_xticks([])
        axs[c].set_yticks([])
    
    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    axs[n//2].set_title(suptitle)
    return fig

from matplotlib import ticker
def plot_from_cooler_and_anchors(coolfile, anc1, anc2, d, axs = None, loopdict=None, cond=None, cols=6, suptitle=None, fgsz=(30*mm, 30*mm)):
    n = 1
    rows = math.ceil(n/cols)
    if axs is None:
        fig, axs = init_subplots_exact(1, 1, fgsz=fgsz, dpi=200, as_list=False)
        axs = np.ravel([axs])
    anc1, anc2 = extend_l(anc1, d), extend_l(anc2, d)
    m = fetch_from_cooler(coolfile, anc1, anc2)

    L, R = anc2[1:]
    T, B = anc1[1:]

    if intersects_diagonal(anc1, anc2):
        axs[0].matshow(m, cmap='gist_heat_r', norm=colors.LogNorm(),
                       extent = [L, R, B, T],
                       )
    else:
        axs[0].matshow(m, cmap='gist_heat_r', 
                       extent = [L, R, B, T])
    if loopdict is not None:
        decorate_ax_with_loops(axs[0], anc1, anc2, loopdict[cond])
    axs[0].set_title(cond)

    ax = axs[0]
    ax.tick_params(axis='both', which='both', right=False, top=False, 
                    labelright=False, labeltop=False, bottom=True, labelbottom=True,
                    )  # Remove tick labels from top and right axes
    ax.ticklabel_format(axis='both', style='sci', scilimits=(6, 6))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"

    ax.xaxis.get_offset_text().set_visible(False)  # Hide the offset text on the x-axis
    ax.yaxis.get_offset_text().set_visible(False)  # Hide the offset text on the y-axis

    chrom = anc1[0]
    ax.set_xlabel(f"Mb, chr{chrom}")  # Set the x-axis label for the last subplot
    ax.set_ylabel(f"Mb, chr{chrom}")  # Set the x-axis label for the last subplot
    ax.xaxis.set_label_coords(0.5, -0.12)  # Adjust the position of the x-axis label

    ax.set_xticks([L, (L+R)//2, R])
    ax.set_yticks([B, (B+T)//2, T])

    # return fig


def get_pval_string_from_pval(pval):
    if pval > .05: 
        pval_string = 'NS'
    else: 
        rest_pval = math.ceil(np.abs((np.log10(pval))))

        front = int(np.round(pval*(10**rest_pval)))
        # pval_string = f'{front}e-{int(rest_pval)}'
        pval_string = f'{pval:.0e}'.lower()
    return pval_string


from aux_functions import grange_to_tuple, make_int
def decorate_ax_with_loops(ax, tuple1, tuple2, loop_bedtool, facecolor='white', alpha=.4, s=200, linewidth=1, **kwargs):
    chrom, s1, e1 = tuple1
    chrom, s2, e2 = tuple2
    reg = [[chrom, s1, e2]]
    for i in loop_bedtool.pair_to_bed(reg, type='ispan'):
        l1 = float(i[1])
        l1e = float(i[1+1])
        l2 = float(i[4])
        l2e = float(i[4+1])
        if ((l1 > s1) and (l1 < e1)) and ((l2 > s2) and (l2 < e2)):
            x = (l1+l1e)//2
            y = (l2+l2e)//2
            ax.scatter(y, x, marker='s', alpha = alpha, c=facecolor, edgecolor='black', s=s, linewidth=linewidth)

import math
import numpy as np
import matplotlib.pyplot as plt
def init_subplots(n, rows, fgsz=(4, 4), space = .2, **kwargs):
    cols = math.ceil(n/rows)
    fig, axs = plt.subplots(rows, cols, figsize=(cols*fgsz[0], rows*fgsz[1]), **kwargs)
    fig.subplots_adjust(hspace=space, wspace=space)
    axs = np.ravel(axs)
    return fig, axs


def init_subplots_exact(n, rows, fgsz=(4, 4), space=1.2, xspace=None, yspace=None, dpi=30, as_list=False, sharex=False, sharey=False, x_ratios=None, y_ratios = None, **kwargs):
    cols = math.ceil(n / rows)
    fig = plt.figure(figsize=fgsz, dpi=dpi)
    axs = []
    first_col_axs = [None] * cols  # First axes in each column for sharey
    first_row_axs = [None] * rows  # First axes in each row for sharex

    for i in range(n):
        row = i // cols
        col = i % cols
        if y_ratios is None:
            yscale = 1
            tot_yscale = 1*row
        else:
            yscale = y_ratios[row]
            tot_yscale = np.sum(y_ratios[:row+1])
        if x_ratios is None:
            xscale = 1
            tot_xscale = 1*col
        else:
            xscale = x_ratios[col]
            tot_xscale = np.sum(x_ratios[:col])
        if xspace is None:
            xspace = space
        if yspace is None:
            yspace = space
        if i == 0:  # The first subplot
            ax = fig.add_axes([tot_xscale * xspace, (n - 1 - tot_yscale) * yspace, xscale, yscale], **kwargs)
            first_col_axs[col] = ax
            first_row_axs[row] = ax
        else:
            ax = fig.add_axes([tot_xscale * xspace, (n - 1 - tot_yscale) * yspace, xscale, yscale], 
                              sharex=first_col_axs[col] if sharex else None, 
                              sharey=first_row_axs[row] if sharey else None, **kwargs)
            if first_col_axs[col] is None:
                first_col_axs[col] = ax
            if first_row_axs[row] is None:
                first_row_axs[row] = ax
        axs.append(ax)

    axs = np.ravel(axs)

    # Remove extra axes if needed
    for ax in axs[n:]:
        ax.remove()

    if (len(axs) == 1) and (not as_list):
        return fig, axs[0]
    else:
        return fig, axs





from matplotlib import transforms

from mpl_toolkits import axes_grid1
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from aux_functions import *

def make_off_diag_plot(to_plot, cooldict, all_ind_to_region, delta=0, vmin = 1e-5, vmax=5e-3, useSigma=True, figsize=(3, 3),
                        yticksToPlot = [0, 1e-3, 2e-3, 3e-3], add_chr=True, key2name={}):
    (ind1, ind2) = to_plot['inds']
    chrom1, s1, _ = all_ind_to_region[ind1[0]]
    chrom1, _, e1 = all_ind_to_region[ind1[1]]
    
    chrom2, s2, _ = all_ind_to_region[ind2[0]]
    chrom2, _, e2 = all_ind_to_region[ind2[1]]
    if add_chr == False:
        vdict = get_merged_mats(cooldict, chrom1, s1, e1, chrom2, s2, e2, delta)
    else:
        vdict = get_merged_mats(cooldict, 'chr' + chrom1, s1, e1, 'chr' + chrom2, s2, e2, delta)
    res = list(cooldict.values())[0].info['bin-size']
    n_conds = len(vdict)
    figlist = []
    mlist = []
    
    fig, axs = init_subplots(n_conds, 1, fgsz = figsize)
    for i, cond in enumerate(vdict.keys()):
        ax = axs[i]
        m = vdict[cond]
        m[np.isnan(m)] = 0
        if useSigma:
            m = scipy.ndimage.gaussian_filter(m, sigma=.5) 
        mlist.append(m)
        pos = ax.matshow(m, cmap=cm.gist_heat_r, vmin=vmin, vmax=vmax)
        n = len(m)
        a = ax

        xticklabels = [s2//1_000_000]
        xticks = [0]
        yticklabels = [s1//1_000_000]
        yticks = [0]
        for ind, lab in zip(to_plot['xlabel'], to_plot['xlabels']):
            x = (ind-ind2[0])*250_000//res+2
            xticks.append(x)
            xticklabels.append(lab)
        for ind, lab in zip(to_plot['ylabel'], to_plot['ylabels']):
            x = (ind-ind1[0])*250_000//res+2
            yticks.append(x)
            yticklabels.append(lab)            
        xticklabels.append(e2//1_000_000)
        xticks.append(len(m))

        yticklabels.append(e1//1_000_000)
        yticks.append(len(m))
        xticks = arr(xticks) - .5
        yticks = arr(yticks) - .5
        
        a.set_xticks(xticks)
        a.set_xticklabels(xticklabels)
        a.set_yticks(yticks)
        a.set_yticklabels(yticklabels)
        a.set_xlabel(f"Mb (chr{chrom1})")
        a.set_ylabel(f"Mb (chr{chrom2})")
        a.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
        a.set_title(f'{key2name.get(cond, cond)}')
        if i == n_conds-1:
            newax = inset_axes(a, width = '4%', height = '100%', 
                        loc = 'lower left', bbox_to_anchor=(1.1, 0, 1, 1), bbox_transform=a.transAxes, borderpad=0)
            plt.colorbar(pos, cax=newax)
            newax.set_yticks(yticksToPlot)
            newax.set_ylim([0, vmax])
        ax.grid(False)
    fig.tight_layout()
    return figlist, mlist



def make_comprehensive_plot(i1, i2, d=40, res=50_000, poiss_vmax=20, extend_by = 250_000):
    grange1 = (chrom, i1*res, i1*res+res)
    grange2 = (chrom, i2*res, i2*res+res)
    grange1, grange2 = map(lambda x: extend_l(x, extend_by), [grange1, grange2])
    collapsemat, filtmat, logp_mat, peak_row_starts, peak_col_starts = get_megaloops_at_grange(grange1, grange2, extend_by = 0)
    L, R = grange2[1:]
    T, B = grange1[1:]
    fig, axs = init_subplots(3, 1,)
    ax = axs[0]
    test = plot_mat(collapsemat, ax=ax, cmap='gist_heat_r', extent = [L, R, B, T], zorder=1);
    ax = axs[1]
    test = plot_mat(filtmat, ax=ax, cmap='gist_heat_r', extent = [L, R, B, T], zorder=1);
    ax = axs[2]
    test = plot_mat(logp_mat, ax=ax, cmap='bwr', extent = [L, R, B, T], zorder=1);
    
    ax = axs[0]
    y, x = peak_row_starts, peak_col_starts
    ax.scatter(x, y, facecolor='none', edgecolor='blue', linewidth=2, s=280, label='Megaloops')
    ax.legend()
    titles = ['Raw Data', 'Smoothed Data', 'Poisson P-value']
#     chrom1, chrom2 = all_ind_to_region_50kb[i1][0], all_ind_to_region_50kb[i2][0]
    for c, ax in enumerate(axs):
        ax.set_title(titles[c])
        ax.tick_params(axis='both', which='both', right=False, top=False, 
                       labelright=False, labeltop=False, bottom=True, labelbottom=True,
                      )  # Remove tick labels from top and right axes
        ax.ticklabel_format(axis='both', style='sci', scilimits=(6, 6))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.0f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.0f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"

        ax.xaxis.get_offset_text().set_visible(False)  # Hide the offset text on the x-axis
        ax.yaxis.get_offset_text().set_visible(False)  # Hide the offset text on the y-axis

        ax.set_xlabel(f"Mb, Chr{chrom}")  # Set the x-axis label for the last subplot
        ax.set_ylabel(f"Mb, Chr{chrom}")  # Set the x-axis label for the last subplot
        ax.xaxis.set_label_coords(0.5, -0.12)  # Adjust the position of the x-axis label

        ax.set_xticks([L, (L+R)//2, R])
        ax.set_yticks([B, (B+T)//2, T])
    plt.tight_layout()
    return collapsemat    


from matplotlib_venn import venn2, venn3
def my_venn2(*args, **kwargs):
    v1 = venn2(*args, **kwargs)
    
    for v in [v1]:
        for text in v.subset_labels:
            if text is None:
                continue
            else:
                text.set_fontsize(6)
        for text in v.set_labels:
            text.set_fontsize(6)
    

def my_venn3(*args, **kwargs):
    v1 = venn3(*args, **kwargs)
    for v in [v1]:
        for text in v.subset_labels:
            if text is None:
                continue
            else:
                text.set_fontsize(6)
        for text in v.set_labels:
            text.set_fontsize(6)
    


def plot_overlaps(bedtool_dict, baseline, jumpsize = 5_000, jumps = np.arange(-10, 10, 1)):
    jumpsizes = jumpsize*jumps
    for name, bedtool in bedtool_dict.items():
        ps = []
        for shift in jumpsizes:
            shifted_bedtool = bedtool.shift(s=shift, genome='mm10')
            p = bedprop(shifted_bedtool, baseline)  
            ps.append(p)
        plt.plot(jumpsizes, ps, label=name)
    plt.legend()
    plt.xlabel('Shift (bp)')
    plt.ylabel('Proportion Overlapping')
    plt.title('Proportion Overlapping vs. Shift')


import matplotlib.gridspec as gridspec
def plot_with_break(X, break_at = 100, ymax = 200, hspace=.2):
    # Set up the gridspec to define the layout of the axes
    fig = plt.figure(figsize=(36*mm, 18*mm), dpi=300)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], hspace=hspace)

    # Create two axes
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    # Plot the same data on both axes
    sns.histplot(X, zorder=3, ax=ax0, stat='count')
    sns.histplot(X, zorder=3, ax=ax1, stat='count')

    # Set the y-axis limits
    ax0.set_ylim(break_at, ymax)
    ax1.set_ylim(0, break_at)
    ax1.set_yticks([0, break_at//2, break_at])
    # Hide the spines between ax and ax2
    ax0.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # ax0.xaxis.tick_top()
    ax0.tick_params(labeltop=False, bottom=False)  # don't put tick labels at the top

    # Set labels and title
    for ax in [ax0, ax1]:
        plt.sca(ax)
        plt.yticks(fontsize=4)

    # Draw diagonal lines on the broken axis
    d = .015  # size of break mark
    kwargs = dict(transform=ax0.transAxes, color='k', clip_on=False)
    ax0.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax0.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    return (fig, (ax0, ax1))




def make_genome_browser_plot(place, peak_set, foxp3_bwdict, histone_bwdict, ylimdict, plot_peaks=True, **kwargs):
    n = len(foxp3_bwdict) + len(histone_bwdict) + 1
    place = place
    fig, axs = plt.subplots(n, 1, figsize=(12, 15), sharex=True)
    ax = axs[0]
    plot_one_GTF_on_axis(place, ax)
    ax.set_yticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(bottom=False)

    for c, name in enumerate(foxp3_bwdict):
        ax = axs[c + 1]
        plot_bigwig_on_axis(place, 
                            foxp3_bwdict, ax, name, ylimdict=ylimdict, label_for_title=name, alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))


    for c, name in enumerate(histone_bwdict):
        ax = axs[c + 1 + len(foxp3_bwdict)]
        plot_bigwig_on_axis(place, 
                            histone_bwdict, ax, name, ylimdict=ylimdict, label_for_title=name,
                            alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))
    if plot_peaks:
        peaks = peak_set.intersect([place], u=True).sort().merge()
        for ax in axs[1:]:
            seen = set()
            for peak in peaks:
                s, e = make_int(peak)[1:3]
                if (s, e) in seen:
                    continue
                else: 
                    seen.add((s, e))
                ax.axvspan(s, e, color='red', alpha=.2, linewidth=0)
    for ax in axs:
        ax.grid(False)    


