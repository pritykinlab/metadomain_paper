from get_focal_contacts_with_peak_prominence import *
from aux_functions import *

def get_megaloops_at_grange_with_collapsing(l1, l2, extend_by = 250_000, logp_co=8,
                                            merged_mat = 'r=0.7_noArmCl13.mcool', resolution=5000,
                                            filter_width=3,
                                            filter_n=35):
    pref = '/Genomics/argo/users/gdolsten/pritlab/mega_tcell_dataset/final_mega_coolfiles/'

    cool = cooler.Cooler(pref + merged_mat + f"::/resolutions/{resolution}")
    l1, l2 = map(add_chr, [l1, l2])
    l1, l2 = map(lambda x: extend_l(x, extend_by), [l1, l2])
    l1, l2 = map(notexpanded_tuple_to_grange, [l1, l2])
    try:
        collapsemat = cool.matrix(balance=False).fetch(grange_to_tuple(l1), grange_to_tuple(l2))
    except Exception as e:
#         print(e.__class__.__name__)
#         if 'Genomic region out of bounds'
        return [None, None, ['Error'], ['Error']]
    w1 = cool.bins().fetch(grange_to_tuple(l1))['weight']
    w2 = cool.bins().fetch(grange_to_tuple(l2))['weight']
    w1[w1 < .0025] = np.nan
    w2[w2 < .0025] = np.nan

    w = np.outer(w1, w2)
    w  = w / np.nanmean(w)
    collapsemat = collapsemat.astype(float)*w
    filtmat, logp_mat, loop_rows, loop_cols, _, _, resultsdict = call_peak_with_poisson_and_peak_cutoff(collapsemat, logp_co=logp_co, filter_n=filter_n, 
                                                                                                        filter_width=filter_width,
                                                                                                        frac_min_valid=.6)
    megaloops_of_interest = (logp_mat > logp_co) & (resultsdict['counts_inner'] > 50) 
    megaloops_of_interest = scipy.ndimage.binary_dilation(megaloops_of_interest)
    outer_edges_mask = np.ones_like(megaloops_of_interest)
    outer_edges_mask[2:-2, 2:-2] = 0
    megaloops_of_interest[outer_edges_mask] = 0

    label_mat, _ = scipy.ndimage.label(megaloops_of_interest)
    obj_locations = scipy.ndimage.find_objects(label_mat)
    collapsed_logp_mat = np.zeros(megaloops_of_interest.shape)
    for c, i in enumerate(obj_locations):
        l, r = i
        # Get OE matrix
        submat = filtmat[l, r].copy() * (label_mat[l, r] == (c+1))
        S, E = np.unravel_index(np.nanargmax(np.ravel(submat)), submat.shape)
        # plt.matshow(submat)
        # plt.title([S, E, l.start, r.start, l, r])
        s1, s2 = l.start, r.start
        collapsed_logp_mat[s1+S, s2+E] = logp_mat[s1+S, s2+E].copy()
    
    s_row = int(grange_to_tuple(l1)[1])
    s_col = int(grange_to_tuple(l2)[1])

    peak_row_inds, peak_col_inds = np.where(collapsed_logp_mat)
    peak_row_starts = s_row + peak_row_inds*resolution
    peak_col_starts = s_col + peak_col_inds*resolution
    return collapsemat, filtmat, logp_mat, peak_row_starts, peak_col_starts

def get_n_foci(reses):
    ls = []
    for i in reses:
        if 'Error' in i[-1]:
            ls.append(np.nan)
        else:
            ls.append(len(i[-1]))
    return ls

def make_foci_count_df(treg_mat, tcon_mat, ls):
    xs, ys = np.where((treg_mat > 0) | (tcon_mat > 0))
    megaloop_pval = treg_mat[xs, ys]
    megaloop_pval_2 = tcon_mat[xs, ys]

    count_df = pd.DataFrame(megaloop_pval, columns=['Treg Pval'])
    count_df['Tcon Pval'] = arr(megaloop_pval_2)
    count_df['Focal peaks'] = arr(ls)
    count_df.fillna(0, inplace=True)
    return count_df


from plotting_functions import init_subplots, plot_mat
import matplotlib.pyplot as plt
from matplotlib import ticker
from plotting_functions import *

def make_comprehensive_plot_grange_collapsing(grange1, grange2, d=40, res=5_000, poiss_vmax=20, extend_by = 250_000, logp_co=8, s = 280,
                                              linewidth=2, fgsz=(30*mm, 30*mm)):
    grange1, grange2 = map(lambda x: extend_l(x, extend_by), [grange1, grange2])
    collapsemat, filtmat, logp_mat, peak_row_starts, peak_col_starts = get_megaloops_at_grange_with_collapsing(grange1, grange2, extend_by = 0, logp_co=logp_co)
    L, R = grange2[1:]
    T, B = grange1[1:]
    fig, axs = init_subplots_exact(3, 1, fgsz=fgsz, dpi = 200, xspace=1.4)
    ax = axs[0]
    test = plot_mat(collapsemat, ax=ax, cmap='gist_heat_r', extent = [L, R, B, T], zorder=1)
    ax = axs[1]
    test = plot_mat(filtmat, ax=ax, cmap='gist_heat_r', extent = [L, R, B, T], zorder=1)
    ax = axs[2]
    test = plot_mat(logp_mat, ax=ax, cmap='bwr', extent = [L, R, B, T], zorder=1, vmin=0, vmax=16)
    
    ax = axs[0]
    y, x = peak_row_starts, peak_col_starts
    ax.scatter(x+.5*res, y+.5*res, facecolor='none', edgecolor='blue', linewidth=linewidth, s=s)
    ax.scatter([], [], facecolor='none', edgecolor='blue', linewidth=1, s=20, label='Foci')
    ax.legend()
    titles = ['Compendium (5kb)', 'Smoothed Data', 'Poisson P-value']

    chrom = grange1[0]
    for c, ax in enumerate(axs):
        ax.set_title(titles[c])
        ax.tick_params(axis='both', which='both', right=False, top=False, 
                       labelright=False, labeltop=False, bottom=True, labelbottom=True,
                      )  # Remove tick labels from top and right axes
        ax.ticklabel_format(axis='both', style='sci', scilimits=(6, 6))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"

        ax.xaxis.get_offset_text().set_visible(False)  # Hide the offset text on the x-axis
        ax.yaxis.get_offset_text().set_visible(False)  # Hide the offset text on the y-axis

        ax.set_xlabel(f"Mb, chr{chrom}")  # Set the x-axis label for the last subplot
        ax.set_ylabel(f"Mb, chr{chrom}")  # Set the x-axis label for the last subplot
        ax.xaxis.set_label_coords(0.5, -0.12)  # Adjust the position of the x-axis label

        ax.set_xticks([L, (L+R)//2, R])
        ax.set_yticks([B, (B+T)//2, T])
    plt.tight_layout()
    return fig, collapsemat

print(1)
def make_simple_plot_grange_collapsing(grange1, grange2, d=40, res=5_000, poiss_vmax=20, extend_by = 250_000, logp_co=8, s = 280,
                                              linewidth=2, fgsz=(40*mm, 40*mm), ignore_set={}):
    grange1, grange2 = map(lambda x: extend_l(x, extend_by), [grange1, grange2])
    collapsemat, filtmat, logp_mat, peak_row_starts, peak_col_starts = get_megaloops_at_grange_with_collapsing(grange1, grange2, extend_by = 0, logp_co=logp_co)
    L, R = grange2[1:]
    T, B = grange1[1:]
    fig, axs = init_subplots_exact(1, 1, fgsz=fgsz, dpi = 200, xspace=1.4, as_list=True)
    ax = axs[0]
    test = plot_mat(collapsemat, ax=ax, cmap='gist_heat_r', extent = [L, R, B, T], zorder=1)

    ax = axs[0]
    y, x = peak_row_starts, peak_col_starts
    ax.scatter(x+.5*res, y+.5*res, facecolor='none', edgecolor='blue', linewidth=linewidth, s=s)
    ax.scatter([], [], facecolor='none', edgecolor='blue', linewidth=1, s=20, label='Foci')
    ax.legend()
    titles = ['Compendium (5kb)', 'Smoothed Data', 'Poisson P-value']

    chrom = grange1[0]

    plt.grid(False)
    tmp1 = add_GTF_to_axis(grange2, plt.gca(), gene_name_fontsize=6, ignore_set = ignore_set,
                    ignore_method='grayout_no_text'
                   )
    tmp1.grid(False)

    tmp2 = add_GTF_to_L_axis(grange1, plt.gca(), gene_name_fontsize=6, ignore_set = ignore_set,
                    ignore_method='grayout_no_text', 
                   )
    tmp2.grid(False)


    for c, ax in enumerate(axs):
        ax.set_title(titles[c])
        ax.tick_params(axis='both', which='both', right=False, top=False, 
                       labelright=False, labeltop=False, bottom=True, labelbottom=True,
                      )  # Remove tick labels from top and right axes
        ax.ticklabel_format(axis='both', style='sci', scilimits=(6, 6))

        ax.xaxis.get_offset_text().set_visible(False)  # Hide the offset text on the x-axis
        ax.yaxis.get_offset_text().set_visible(False)  # Hide the offset text on the y-axis

        # ax.set_xlabel(f"Mb, chr{chrom}")  # Set the x-axis label for the last subplot
        ax.xaxis.set_label_coords(0.5, -0.12)  # Adjust the position of the x-axis label

        tmp1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"
        tmp2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x / 1e6)))  # Change x-axis label format to "Mb"

        ax.set_xticks([])
        ax.set_yticks([])
        tmp1.set_xticks([L, (L+R)//2, R])
        tmp2.set_yticks([B, (B+T)//2, T])
        tmp2.set_ylabel(f"Mb, chr{chrom}")  # Set the x-axis label for the last subplot
    plt.tight_layout()
    return fig, [axs, tmp1, tmp2], collapsemat