import scipy
import scipy.stats
import itertools
import matplotlib.pyplot as plt
import numpy as np

from aux_functions import format_pval_as_asterisks

def add_stat_annotation_boxplot_with_hue(ax, data, xcol, ycol, hue, order, hue_order, box_pairs, ymax=.32, delta = .15,
                                         yoff_method = 'add'):
    """ Add statistical annotations for comparing hue within each x."""
    unique_x = data[xcol].unique()
    hue_levels = data[hue].unique()
    
    # Prepare a mapping from hue and x to x-tick positions
    hue_offsets = {hue: (i-1)*.25 for i, hue in enumerate(hue_order)}
    x_positions = {v: k for k, v in enumerate(order)}

    seen = set()
    yoff = 0
    for pair in box_pairs:
        if yoff_method == 'add':
            if pair[0][0] not in seen:
                seen.add(pair[0][0])
                yoff = 0 
            else:
                yoff += delta
        else:
            yoff += delta
        
        print(yoff)
        # Filter data for each pair

        data1 = data[(data[xcol] == pair[0][0]) & (data[hue] == pair[0][1])]
        data2 = data[(data[xcol] == pair[1][0]) & (data[hue] == pair[1][1])]

        # Perform rank sum test
        stat, p_value = scipy.stats.ranksums(data1[ycol], data2[ycol])
    
        # Find positions for annotations
        x1 = x_positions[pair[0][0]] + hue_offsets[pair[0][1]] + .25/2
        x2 = x_positions[pair[1][0]] + hue_offsets[pair[1][1]] + .25/2
        y, h, col = ymax + 0.02, 0.01, 'k'

        y = y * (1 + yoff)
        # Draw the lines and annotations
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.2, c=col)
        ax.text((x1 + x2) * .5, y + h, f'{format_pval_as_asterisks(p_value)}', ha='center', va='bottom', color=col,
               fontsize=4)


def add_stat_annotation_boxplot_no_hue(ax, data, xcol, ycol, order, box_pairs, ymax=.32, delta = .15, h = 0.01,
                                       log=False, x_skip = 1):
    """ Add statistical annotations for comparing hue within each x."""
    unique_x = data[xcol].unique()
    # Prepare a mapping from hue and x to x-tick positions
    x_positions = {v: k*x_skip for k, v in enumerate(order)}

    seen = set()
    yoff = 0
    for pair in box_pairs:

        yoff += delta
        # Filter data for each pair

        data1 = data[(data[xcol] == pair[0])].dropna()
        data2 = data[(data[xcol] == pair[1])].dropna()

        # Perform rank sum test
        stat, p_value = scipy.stats.ranksums(data1[ycol], data2[ycol])
        if log:
            print(f'Comparison of {pair[0]} and {pair[1]}')
            print(f'p-value: {p_value}')
            print(f'stat: {stat}')
        # Find positions for annotations
        x1 = x_positions[pair[0]]
        x2 = x_positions[pair[1]]
        y, h, col = ymax + 0.02, h, 'k'

        y = y * (1 + yoff)
        # Draw the lines and annotations
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=.2, c=col)
        ax.text((x1 + x2) * .5, y + h, f'{format_pval_as_asterisks(p_value)}', ha='center', va='bottom', color=col,
               fontsize=4)
