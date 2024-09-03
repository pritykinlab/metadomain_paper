from plotting_functions	import *
import matplotlib.patheffects as pe
from adjustText import adjust_text

def volcano_plot(lfcs, pvals, names, ax, max_y = 12.5, pco = .01, lfc_co = .25, c = None, vmin=0, vmax=1,
                    numerator = '', denominator = '', title = '', xlabel = '', ylabel = '', rasterized = True, s=4,
                    fontsize = 8, label_pval_cutoff = 3, label_lfc_cutoff = 0, up_color = 'salmon', down_color = 'skyblue',
                    xlim = [-2, 2], ylim = [0, 12.5], ignore_set = set(), alpha=1, do_adjust_text=True,
                    genes_to_label = [], pes=[pe.withStroke(linewidth = .25, foreground="black")]):
    texts = []
    # x, y = joris_lfc_df.loc[names, col], -np.log10(joris_pval_df.loc[names, col]).copy()
    x, y = lfcs, -np.log10(pvals)
    y[y > max_y] = max_y
    colors = np.asarray(['lightgray']*len(x))
    pco = -np.log10(pco)
    up_idx = (y > pco) & (x > lfc_co)
    down_idx = (y > pco) & (x < -lfc_co)
    if c is None:
        colors[up_idx] = up_color
        colors[down_idx] = down_color
    else:
        colors = c
    plt.scatter(x, y, zorder = 3, c = colors, linewidth=0, s = s, rasterized=rasterized, vmin=vmin, vmax=vmax, alpha=alpha)
    plt.ylabel('-$\log_{10}$(FDR)')

    plt.text(.95, .05, s = f'{up_idx.sum()}', transform = ax.transAxes, color=up_color, ha='right', va='bottom', fontsize=fontsize)
    plt.text(0.05, .05, s = f'{down_idx.sum()}', transform = ax.transAxes, color=down_color, ha='left', va='bottom', fontsize=fontsize)

    plt.text(.95, -.2, s = f'{numerator}', transform = ax.transAxes,  ha='right', va='bottom', fontsize=fontsize)
    plt.text(0.05, -.2, s = f'{denominator}', transform = ax.transAxes,  ha='left', va='bottom', fontsize=fontsize)
    texts = []
    for i in np.where((y > label_pval_cutoff) & (x.abs() > label_lfc_cutoff) | (names.isin(genes_to_label)))[0]:
        if x.iloc[i] > 0:
            color = up_color
            textcolor = 'red'
        else:
            color = down_color
            textcolor = 'blue'            
        gene_name = names.iloc[i]
        if (('Gm' in gene_name) or ('Rik' in gene_name) or ('ENSMUS' in gene_name) or (gene_name in ignore_set) or ('-ps' in gene_name)
                or ('None' in gene_name)):
            print(f'ignoring {gene_name}')
            continue
        if np.isnan(x.iloc[i]) or np.isnan(y.iloc[i]):
            continue
        plt.scatter(x.iloc[i], y.iloc[i], zorder = 3.6, color = 'None', edgecolor = 'black', linewidth = .2, s = s)
        texts.append(plt.text(x.iloc[i], y.iloc[i], gene_name, fontsize=fontsize, color = textcolor, 
                             path_effects=pes, zorder = 4,
                              
                             ))
    plt.xlim(xlim)
    plt.ylim(ylim)
    if do_adjust_text:
        adjust_text(texts, ax=plt.gca(), arrowprops=dict(arrowstyle="-", color='black', lw=0.15, alpha=.8), zorder=3.5,
                # explode_radius = 150
                )

import matplotlib.patheffects as pathe
def label_points(x, y, indices, labels, fontsize=6, textcolor='black', s=4, ha='right', linewidth = .2, bold=True):
    texts = []
    for c, i in enumerate(indices):
        gene_name = labels[c]
        plt.scatter(x[i], y[i], zorder = 3.6, color = 'none', edgecolor = 'black', linewidth = linewidth, s = s)
        if bold:
            pe=[pathe.withStroke(linewidth = .25, foreground="black")]
        else:
            pe = []
        texts.append(plt.text(x[i], y[i], gene_name, fontsize=fontsize, color = textcolor, 
                             path_effects=pe, zorder = 4,
                             ha='left',
                             ma = 'left',
                              
                             ))
    adjust_text(texts, ax=plt.gca(), arrowprops=dict(arrowstyle="-", color='black', lw=0.15, alpha=.8), zorder=3.5)