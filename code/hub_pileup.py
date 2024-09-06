import time
import os
from aux_functions import *
from get_focal_contacts_with_peak_prominence import *

def get_inter_hub_peak_hichip(all_inds, all_ind_to_region, merged_clustdict, goodinds, my_treg_comp, merged_inds_to_subset,
                            chrom_to_start, chrom_to_end):
    rows = []
    for ind in all_inds:
        t1 = time.time()
        if not os.path.exists(f'for_snipping/hichip_pileups/ind={ind}_inds_to_compare_with.pkl'):
            continue
        pileup_treg_hichip_lol = load_pickle(f'for_snipping/hichip_pileups/ind={ind}_tTreg_pileup_dict.pkl')
        pileup_cd4sp_hichip_lol = load_pickle(f'for_snipping/hichip_pileups/ind={ind}_CD4SP_pileup_dict.pkl')
        pileup_prec_hichip_lol = load_pickle(f'for_snipping/hichip_pileups/ind={ind}_Precursor_pileup_dict.pkl')
        pileup_inds = load_pickle(f'for_snipping/hichip_pileups/ind={ind}_inds_to_compare_with.pkl')

        pileup_dict = {
            'tTreg ÷ tTreg Prec' : (pileup_treg_hichip_lol, pileup_prec_hichip_lol),
            'tTreg ÷ CD4SP' : (pileup_treg_hichip_lol, pileup_cd4sp_hichip_lol),
            'tTreg Prec ÷ CD4SP' : (pileup_prec_hichip_lol, pileup_cd4sp_hichip_lol),
        }
        for name, (pileup1, pileup2) in pileup_dict.items():
            for c, u in enumerate(merged_inds_to_subset):
                chrom = all_ind_to_region[ind][0]
                s = chrom_to_start[chrom]
                e = chrom_to_end[chrom]
                good = goodinds[(merged_clustdict['all'] == u) & ((goodinds < s) | (goodinds > e))]
                indsoi = np.isin(pileup_inds, good)
                    
                treg_up, tcon_up, _, _, _ = compute_interchromosomal_hub_pileup(pileup1, pileup2, ind, indsoi)
                rows.append([name, treg_up, tcon_up, u, ind])
        
            generic_a_compartment = np.where(my_treg_comp > .9)[0]
            generic_a_compartment = generic_a_compartment[~np.isin(generic_a_compartment, goodinds)]
            generic_a_compartment = generic_a_compartment[(generic_a_compartment < s) | (generic_a_compartment > e)]
            generic_a_compartment = np.isin(pileup_inds, generic_a_compartment)
            treg_up, tcon_up, _, _, _ = compute_interchromosomal_hub_pileup(pileup1, pileup2, ind, generic_a_compartment)
            rows.append([name, treg_up, tcon_up, 'A Compartment', ind])

        loc = ind
        t1f = time.time()
        if loc % 10 == 0:
            print(f"Done with {loc}: {t1f-t1} \n")
    print("Done")
    data = pd.DataFrame(rows, columns=['name', 'treg_up', 'tcon_up', 'hub', 'ind'])
    return data


def get_inds_for_oe_df(all_intra_metadomains, bedtool_dict, all_ind_to_region, chrom_to_start, chrom_to_end, CHROMS_TO_USE):
    l = 0
    index_dict = {}
    for key, bedtool in bedtool_dict.items():
        indices = []
        for chrom in CHROMS_TO_USE:
            s, e = chrom_to_start[chrom], chrom_to_end[chrom]
            sl1 = slice(s, e)
        
            regions = add_chr_to_bedtool(all_ind_to_region[s:e])
            has_bedtool = get_col(regions.intersect(bedtool, c=True), -1).astype(float)
        
            stat5_outer = np.add.outer(has_bedtool, has_bedtool)
        
            idx = (all_intra_metadomains[sl1, sl1] > 0) & (stat5_outer > 0)
            metadomains = np.where(np.triu(idx))
            ind = chrom + "_" + pd.Series((metadomains[0]+s).astype('str')) + "_" + pd.Series((metadomains[1]+s).astype('str'))
            indices.extend(ind)
        index_dict[key] = indices
    return index_dict

def get_treg_inds_for_oe_df(all_intra_metadomains, deseq_effect_mat,
                        chrom_to_start, chrom_to_end, CHROMS_TO_USE, cutoff=4, ):
    indices = []
    for chrom in CHROMS_TO_USE:
        s, e = chrom_to_start[chrom], chrom_to_end[chrom]
        sl1 = slice(s, e)

        metadomains = np.where(np.triu((deseq_effect_mat[sl1, sl1] > cutoff) & (all_intra_metadomains[sl1, sl1] > 0)))
        ind = chrom + "_" + pd.Series((metadomains[0]+s).astype('str')) + "_" + pd.Series((metadomains[1]+s).astype('str'))
        indices.extend(ind)
    return indices




from collections import defaultdict
from make_figure4 import *
import itertools

def pileup_bins_with_hub(df, metadomain_pileup_cooldict, chrom_to_start, inter_and_intra_connections_treg, 
                            inter_and_intra_connections_tcon, resolution_in=250_000, resolution_out = 50_000, 
                            intra=False, inter=True, padding_size=30, fetch_oe=True, skip_metadomains=False):
    all_mat_dict = defaultdict(list)
    all_metadata = []
    chromsoi = sorted(df['chrom'].unique())
    stride = resolution_in // resolution_out
    assert resolution_in % resolution_out == 0

    for chrom1, chrom2 in itertools.product(chromsoi, chromsoi):
        if (chrom1 < chrom2):
            continue
        elif (chrom1 == chrom2) and (intra == False):
            continue
        elif (chrom1 != chrom2) and (inter == False):
            continue
        else:
            idx1 = df['chrom'] == chrom1
            idx2 = df['chrom'] == chrom2
            subdf_X = df[idx1].copy()
            subdf_Y = df[idx2].copy()
            subdf_X['norm_ind'] = subdf_X['ind'] - chrom_to_start[chrom1]
            subdf_Y['norm_ind'] = subdf_Y['ind'] - chrom_to_start[chrom2]
            rows, cols = list(zip(*itertools.product(subdf_X.index, subdf_Y.index)))
            
            for name, cool in metadomain_pileup_cooldict.items():
                assert cool.info['bin-size'] == resolution_out, (cool.info['bin-size'], resolution_out)
                append_mats(subdf_X, subdf_Y, cool, name, rows, cols, all_mat_dict, 
                            padding_size=padding_size, stride=stride, fetch_oe=fetch_oe)

            for row, col in zip(rows, cols):
                row_X = subdf_X.loc[row]
                row_Y = subdf_Y.loc[col]

                row, col = stride*row_X['norm_ind']+padding_size, stride*row_Y['norm_ind']+padding_size
                ind_X, ind_Y = row_X['ind'], row_Y['ind']
                cluster_X, cluster_Y = row_X['cluster'], row_Y['cluster']

                if skip_metadomains:
                    mega_loops_treg = None
                    mega_loops_tcon = None
                else:
                    mega_loops_treg = inter_and_intra_connections_treg[ind_X, ind_Y]
                    mega_loops_tcon = inter_and_intra_connections_tcon[ind_X, ind_Y]

                all_metadata.append([mega_loops_treg, mega_loops_tcon, ind_X, ind_Y, cluster_X, cluster_Y])
        print("Done with", chrom1)
    all_metadata = pd.DataFrame(all_metadata, columns = ['treg_mega', 'tcon_mega', 'ind1', 'ind2', 'cluster1', 'cluster2'])
    for key, v in all_mat_dict.items():
        all_mat_dict[key] = np.stack(all_mat_dict[key])
    return all_mat_dict, all_metadata


def pileup_bins_at_diag(df, metadomain_pileup_cooldict, chrom_to_start, inter_and_intra_connections_treg, 
                            inter_and_intra_connections_tcon, resolution_in=250_000, resolution_out = 50_000, 
                            intra=False, inter=True, padding_size=30, fetch_oe=True, skip_metadomains=False):
    all_mat_dict = defaultdict(list)
    all_metadata = []
    chromsoi = sorted(df['chrom'].unique())
    stride = resolution_in // resolution_out
    assert resolution_in % resolution_out == 0

    for chrom1, chrom2 in itertools.product(chromsoi, chromsoi):
        if (chrom1 < chrom2):
            continue
        elif (chrom1 == chrom2) and (intra == False):
            continue
        elif (chrom1 != chrom2) and (inter == False):
            continue
        else:
            idx1 = df['chrom'] == chrom1
            idx2 = df['chrom'] == chrom2
            subdf_X = df[idx1].copy()
            subdf_Y = df[idx2].copy()
            subdf_X['norm_ind'] = subdf_X['ind'] - chrom_to_start[chrom1]
            subdf_Y['norm_ind'] = subdf_Y['ind'] - chrom_to_start[chrom2]
            rows, cols = subdf_X.index, subdf_Y.index
            
            for name, cool in metadomain_pileup_cooldict.items():
                assert cool.info['bin-size'] == resolution_out, (cool.info['bin-size'], resolution_out)
                append_mats(subdf_X, subdf_Y, cool, name, rows, cols, all_mat_dict, 
                            padding_size=padding_size, stride=stride, fetch_oe=fetch_oe)

            for row, col in zip(rows, cols):
                row_X = subdf_X.loc[row]
                row_Y = subdf_Y.loc[col]

                row, col = stride*row_X['norm_ind']+padding_size, stride*row_Y['norm_ind']+padding_size
                ind_X, ind_Y = row_X['ind'], row_Y['ind']
                cluster_X, cluster_Y = row_X['cluster'], row_Y['cluster']

                if skip_metadomains:
                    mega_loops_treg = None
                    mega_loops_tcon = None
                else:
                    mega_loops_treg = inter_and_intra_connections_treg[ind_X, ind_Y]
                    mega_loops_tcon = inter_and_intra_connections_tcon[ind_X, ind_Y]

                all_metadata.append([mega_loops_treg, mega_loops_tcon, ind_X, ind_Y, cluster_X, cluster_Y])
        print("Done with", chrom1)
    all_metadata = pd.DataFrame(all_metadata, columns = ['treg_mega', 'tcon_mega', 'ind1', 'ind2', 'cluster1', 'cluster2'])
    for key, v in all_mat_dict.items():
        all_mat_dict[key] = np.stack(all_mat_dict[key])
    return all_mat_dict, all_metadata








def pileup_bin_pairs(df1, df2, metadomain_pileup_cooldict, chrom_to_start, inter_and_intra_connections_treg, 
                            inter_and_intra_connections_tcon, resolution_in=250_000, resolution_out = 50_000, 
                            intra=False, inter=True, padding_size=30, fetch_oe=True, skip=1, skip_metadomains=False):
    all_mat_dict = defaultdict(list)
    all_metadata = []
    chromsoi = list(set(list(df1['chrom'].unique()) + list(df2['chrom'].unique())))
    stride = resolution_in // resolution_out
    assert resolution_in % resolution_out == 0

    for chrom1, chrom2 in itertools.product(chromsoi, chromsoi):
        if (chrom1 == chrom2) and (intra == False):
            continue
        elif (chrom1 != chrom2) and (inter == False):
            continue
        else:
            idx1 = df1['chrom'] == chrom1
            idx2 = df2['chrom'] == chrom2
            if (idx1.sum() == 0) or (idx2.sum() == 0):
                continue
            subdf_X = df1[idx1].copy()
            subdf_Y = df2[idx2].copy()
            subdf_X['norm_ind'] = subdf_X['ind'] - chrom_to_start[chrom1]
            subdf_Y['norm_ind'] = subdf_Y['ind'] - chrom_to_start[chrom2]
            rows, cols = list(zip(*itertools.product(subdf_X.index, subdf_Y.index)))
            rows, cols = rows[::skip], cols[::skip]
            for name, cool in metadomain_pileup_cooldict.items():
                assert cool.info['bin-size'] == resolution_out
                append_mats(subdf_X, subdf_Y, cool, name, rows, cols, all_mat_dict, 
                            padding_size=padding_size, stride=stride, fetch_oe=fetch_oe)

            for row, col in zip(rows, cols):
                row_X = subdf_X.loc[row]
                row_Y = subdf_Y.loc[col]

                row, col = stride*row_X['norm_ind']+padding_size, stride*row_Y['norm_ind']+padding_size
                ind_X, ind_Y = row_X['ind'], row_Y['ind']
                cluster_X, cluster_Y = row_X['cluster'], row_Y['cluster']

                if skip_metadomains:
                    mega_loops_treg = None
                    mega_loops_tcon = None
                else:
                    mega_loops_treg = inter_and_intra_connections_treg[ind_X, ind_Y]
                    mega_loops_tcon = inter_and_intra_connections_tcon[ind_X, ind_Y]


                all_metadata.append([mega_loops_treg, mega_loops_tcon, ind_X, ind_Y, cluster_X, cluster_Y])
        print("Done with", chrom1)
    all_metadata = pd.DataFrame(all_metadata, columns = ['treg_mega', 'tcon_mega', 'ind1', 'ind2', 'cluster1', 'cluster2'])
    for key, v in all_mat_dict.items():
        all_mat_dict[key] = np.stack(all_mat_dict[key])
    return all_mat_dict, all_metadata


def pileup_direct_bin_pairs(df1, df2, metadomain_pileup_cooldict, chrom_to_start, inter_and_intra_connections_treg, 
                            inter_and_intra_connections_tcon, resolution_in=250_000, resolution_out = 50_000, 
                            intra=False, inter=True, padding_size=30, fetch_oe=True, skip=1, skip_metadomains=False, log=True):
    all_mat_dict = defaultdict(list)
    all_metadata = []
    chromsoi = list(set(list(df1['chrom'].unique()) + list(df2['chrom'].unique())))
    stride = resolution_in // resolution_out
    assert resolution_in % resolution_out == 0

    for chrom1, chrom2 in itertools.product(chromsoi, chromsoi):
        if (chrom1 == chrom2) and (intra == False):
            continue
        elif (chrom1 != chrom2) and (inter == False):
            continue
        else:
            idx1 = df1['chrom'] == chrom1
            idx2 = df2['chrom'] == chrom2
            idx = idx1 & idx2
            if (idx1.sum() == 0) or (idx2.sum() == 0):
                continue
            subdf_X = df1[idx].copy()
            subdf_Y = df2[idx].copy()
            subdf_X['norm_ind'] = subdf_X['ind'] - chrom_to_start[chrom1]
            subdf_Y['norm_ind'] = subdf_Y['ind'] - chrom_to_start[chrom2]
            rows, cols = subdf_X.index, subdf_Y.index

            rows, cols = rows[::skip], cols[::skip]
            for name, cool in metadomain_pileup_cooldict.items():
                assert cool.info['bin-size'] == resolution_out
                append_mats(subdf_X, subdf_Y, cool, name, rows, cols, all_mat_dict, 
                            padding_size=padding_size, stride=stride, fetch_oe=fetch_oe, log=log)

            for row, col in zip(rows, cols):
                row_X = subdf_X.loc[row]
                row_Y = subdf_Y.loc[col]

                row, col = stride*row_X['norm_ind']+padding_size, stride*row_Y['norm_ind']+padding_size
                ind_X, ind_Y = row_X['ind'], row_Y['ind']
                cluster_X, cluster_Y = row_X['cluster'], row_Y['cluster']

                if skip_metadomains:
                    mega_loops_treg = None
                    mega_loops_tcon = None
                else:
                    mega_loops_treg = inter_and_intra_connections_treg[ind_X, ind_Y]
                    mega_loops_tcon = inter_and_intra_connections_tcon[ind_X, ind_Y]


                all_metadata.append([mega_loops_treg, mega_loops_tcon, ind_X, ind_Y, cluster_X, cluster_Y])
        print("Done with", chrom1)
    all_metadata = pd.DataFrame(all_metadata, columns = ['treg_mega', 'tcon_mega', 'ind1', 'ind2', 'cluster1', 'cluster2'])
    for key, v in all_mat_dict.items():
        all_mat_dict[key] = np.stack(all_mat_dict[key])
    return all_mat_dict, all_metadata



def append_mats_OG(subdf_X, subdf_Y, cool, name, rows, cols, mat_dict, padding_size = 30, fetch_oe=True, pc=1e-4, stride=5,):
    chrom1, chrom2 = subdf_X['chrom'].iloc[0], subdf_Y['chrom'].iloc[0]
    padded_mat = fetch_mat(cool, chrom1, chrom2, padding_size=padding_size, fetch_oe=fetch_oe)
    for row, col in zip(rows, cols):
        row_X = subdf_X.at[row, 'norm_ind']
        row_Y = subdf_Y.at[col, 'norm_ind']
        row, col = stride*row_X+padding_size, stride*row_Y+padding_size
        submat = padded_mat[row-padding_size:row+padding_size+1, 
                            col-padding_size:col+padding_size+1]
        assert submat.shape == (2*padding_size+1, 2*padding_size+1), (submat.shape, row, col, padded_mat.shape, row-padding_size, row+padding_size, padding_size)
        mat_dict[name].append(submat)


from construct_oe_mat import fetch_mat
def append_mats(
        subdf_X, subdf_Y, cool, name, rows, cols, mat_dict, padding_size = 30, fetch_oe=True, pc=1e-4, stride=5, log=True,
                ):

    chrom1, chrom2 = subdf_X['chrom'].iloc[0], subdf_Y['chrom'].iloc[0]
    padding_to_center = stride // 2
    padded_mat = fetch_mat(cool, chrom1, chrom2, padding_size=padding_size, fetch_oe=fetch_oe, pc=pc, log=log)
    
    row_Xs = subdf_X.loc[list(rows), 'norm_ind']*stride + padding_size
    row_Ys = subdf_Y.loc[list(cols), 'norm_ind']*stride + padding_size
    for row, col in zip(row_Xs, row_Ys):
        submat = padded_mat[row-padding_size:row+padding_size+1, 
                            col-padding_size:col+padding_size+1][padding_to_center:,
                                                                 padding_to_center:]
        assert submat.shape == (2*padding_size+1 - padding_to_center, 2*padding_size+1 - padding_to_center), (submat.shape, row, col, padded_mat.shape, row-padding_size, row+padding_size, padding_size)
        mat_dict[name].append(submat)


# def append_mats_fast(
#         subdf_X, subdf_Y, cool, name, rows, cols, mat_dict, padding_size = 30, fetch_oe=True, pc=1e-4, stride=5,
#                 ):

#     chrom1, chrom2 = subdf_X['chrom'].iloc[0], subdf_Y['chrom'].iloc[0]
#     padded_mat = fetch_mat(cool, chrom1, chrom2, padding_size=padding_size, fetch_oe=fetch_oe)
#     row_Xs = subdf_X.loc[list(rows), 'norm_ind']*stride+padding_size
#     row_Ys = subdf_Y.loc[list(cols), 'norm_ind']*stride+padding_size

#     row_starts = row_Xs - padding_size
#     col_starts = row_Ys - padding_size

#     # Use broadcasting to generate grids of indices
#     all_rows = row_starts[:, None] + np.arange(2 * padding_size + 1)
#     all_cols = col_starts[:, None] + np.arange(2 * padding_size + 1)

#     # Extract the windows
#     windows = padded_mat[all_rows[:, :, None], all_cols[:, None, :]]
#     return windows

def plot_cluster_pileups_from_result(key, mat_dict, metadata, clusters, res=50_000, columns_to_names={}, row_colors_dict={},
                                     vmin = .5, vmax = 1, s1=0, s2=-1, center_method='center_only', center=2,
                                     show_filts=False, fgsz=(20*mm, 20*mm), log=False, cross_plot=True, xspace=1.2,
                                     cliplo=-1, cliphigh=10, delta_co=0):
    figs = []
    mats = mat_dict[key]
    all_results = {}
    us = np.unique(clusters)
    n = len(us)
    fig, axs = init_subplots_exact(n, 1, fgsz=fgsz, dpi = 100, as_list=True, xspace=xspace)
    for c, (u1) in enumerate(us):

        idx = (metadata['cluster1'] == u1) & (metadata['cluster2'] == u1)
        if isinstance(vmin, list):
            vmin_ = vmin[c]
        else:
            vmin_ = vmin
        if isinstance(vmax, list):
            vmax_ = vmax[c]
        else:   
            vmax_ = vmax

        mat, colorbar, results = plot_pileup_mat(mats[idx], axs[c], vmax=vmax_, vmin=vmin_, res=res, s1=s1, s2=s2,
                                                 method=center_method, center=center, show_filts=show_filts, log=log,
                                                 cliplo=cliplo, cliphigh=cliphigh, delta_co=delta_co)
        plt.title(f'{columns_to_names.get(u1, u1).replace(newline, " ")}', color = row_colors_dict.get(columns_to_names.get(u1, u1),
                                                                                                       'black')
                  )
        all_results[u1] = results
    cbar_ax = axs[-1].inset_axes(transform = axs[c].transAxes, bounds = (1.1, 0, .05, 1))
    plt.colorbar(colorbar, cax=cbar_ax, shrink=.5)
    cbar_ax.grid(False)
    fig.add_axes(cbar_ax)
    plt.sca(cbar_ax)
    plt.yticks(fontsize=4)

    mb = mat.shape[1]*res // 2 / 1e6
    for c, a in enumerate(axs):
        cutoff = np.round(mb, 1)
        a.grid(False)
        a.set_xticks([-mb, 0, mb])
        if c == 0:
            a.set_yticklabels([-cutoff, "Anchor1", cutoff])
        else:
            a.set_yticklabels([])
        a.set_yticks([-mb, 0, mb])
        a.set_xticklabels([-cutoff, "Anchor2", cutoff])
        a.tick_params(labeltop = False, top = False, labelbottom = True, bottom = True)
        a.get_yticklabels()[1].set_rotation(90)
        a.get_yticklabels()[1].set_va('center')
    figs.append(fig)
    mats = mat_dict[key]
    if cross_plot:
        fig, axs = init_subplots_exact(n, 1, fgsz=(20*mm, 20*mm), dpi = 100, as_list=True)
        for c, (u1, u2) in enumerate(itertools.combinations(clusters, 2)):
            idx1 = (metadata['cluster1'] == u1) & (metadata['cluster2'] == u2)
            idx2 = (metadata['cluster1'] == u2) & (metadata['cluster2'] == u1)
            idx = idx1 | idx2
            if isinstance(vmin, list):
                vmin_ = vmin[c+3]
            else:
                vmin_ = vmin
            if isinstance(vmax, list):
                vmax_ = vmax[c+3]
            else:   
                vmax_ = vmax

            mat, colorbar, _ = plot_pileup_mat(mats[idx], axs[c], vmax=vmax_, vmin=vmin_, res=res, s1=s1, s2=s2,
                                            method=center_method, center=center, delta_co=delta_co)
            name1 = columns_to_names.get(u1, u1).replace(newline, "").replace("SE", "").replace("(", "").replace(")", "")
            name2 = columns_to_names.get(u2, u2).replace(newline, "").replace("SE", "").replace("(", "").replace(")", "")
            plt.title(f'{name1} to {name2}')
        
        mb = mat.shape[1]*res //2 / 1e6
        for c, a in enumerate(axs):
            cutoff = np.round(mb, 1)
            a.grid(False)
            a.set_xticks([-mb, 0, mb])
            if c == 0:
                # a.set_yticklabels([-cutoff, name, cutoff])
                a.set_yticklabels([-cutoff, "Anchor1", cutoff])
            else:
                a.set_yticklabels([])
            a.set_yticks([-mb, 0, mb])
            a.set_xticklabels([-cutoff, "Anchor2", cutoff])
            a.tick_params(labeltop = False, top = False, labelbottom = True, bottom = True)
            a.get_yticklabels()[1].set_rotation(90)
            a.get_yticklabels()[1].set_va('center')
        
        cbar_ax = axs[-1].inset_axes(transform = axs[c].transAxes, bounds = (1.1, 0, .05, 1))
        plt.colorbar(colorbar, cax=cbar_ax, shrink=.5)
        cbar_ax.grid(False)
        fig.add_axes(cbar_ax)
        plt.sca(cbar_ax)
        plt.yticks(fontsize=4)

        figs.append(fig)
    return figs, all_results

def plot_ind_pileups_from_result(key, mat_dict, metadata, clusters, ind, res=50_000, columns_to_names={}, row_colors_dict={},
                                     vmin = .5, vmax = 1, s1=0, s2=-1, center_method='None', center=2,
                                     show_filts=False, delta_co=0):
    figs = []
    mats = mat_dict[key]
    all_results = {}
    fig, axs = init_subplots_exact(3, 1, fgsz=(20*mm, 20*mm), dpi = 100)
    for c, (u1) in enumerate(np.unique(clusters)):
        idx1 = (metadata['ind1'] == ind) & (metadata['cluster2'] == u1)
        idx2 = (metadata['ind2'] == ind) & (metadata['cluster1'] == u1).values

        L_mat = mats[idx1]
        R_mat = np.rot90(mats[idx2], k=3, axes=(1, 2))
        # Rotate R_mat
        concat_mats = np.concatenate([L_mat, R_mat], axis=0)
        if isinstance(vmin, list):
            vmin_ = vmin[c]
        else:
            vmin_ = vmin
        if isinstance(vmax, list):
            vmax_ = vmax[c]
        else:   
            vmax_ = vmax
        mat, colorbar, results = plot_pileup_mat(concat_mats, axs[c], vmax=vmax_, vmin=vmin_, res=res, s1=s1, s2=s2,
                                                 method=center_method, center=center, show_filts=show_filts,
                                                 delta_co=delta_co)    
        plt.title(f'{columns_to_names.get(u1, u1).replace(newline, " ")}', color = row_colors_dict[columns_to_names[u1]])
        all_results[u1] = results
    cbar_ax = axs[-1].inset_axes(transform = axs[c].transAxes, bounds = (1.1, 0, .05, 1))
    plt.colorbar(colorbar, cax=cbar_ax, shrink=.5)
    cbar_ax.grid(False)
    fig.add_axes(cbar_ax)
    plt.sca(cbar_ax)
    plt.yticks(fontsize=4)

    mb = mat.shape[1]*res // 2 / 1e6
    for c, a in enumerate(axs):
        cutoff = np.round(mb, 1)
        a.grid(False)
        a.set_xticks([-mb, 0, mb])
        if c == 0:
            a.set_yticklabels([-cutoff, "Anchor1", cutoff])
        else:
            a.set_yticklabels([])
        a.set_yticks([-mb, 0, mb])
        a.set_xticklabels([-cutoff, "Anchor2", cutoff])
        a.tick_params(labeltop = False, top = False, labelbottom = True, bottom = True)
        a.get_yticklabels()[1].set_rotation(90)
        a.get_yticklabels()[1].set_va('center')

    figs.append(fig)
    return figs, all_results


def plot_all_pairs(key, mat_dict, metadata, clusters, res=50_000, columns_to_names={}, row_colors_dict={},
                                     vmin = .5, vmax = 1, s1=0, s2=-1, center_method='None', center=2,
                                     show_filts=False, delta_co=0):
    figs = []
    mats = mat_dict[key]
    all_results = {}
    n = len(np.unique(clusters))
    fig, axs = init_subplots_exact(n*n, n, fgsz=(20*mm, 20*mm), dpi = 100)
    print(len(axs))
    for c, (u1, u2) in enumerate(itertools.product(np.unique(clusters), np.unique(clusters))):
        idx1 = (metadata['cluster1'] == u1) & (metadata['cluster2'] == u2) & ((metadata['ind1'] - metadata['ind2']).abs() > 10)
        idx2 = (metadata['cluster2'] == u1) & (metadata['cluster1'] == u2) & ((metadata['ind1'] - metadata['ind2']).abs() > 10)

        mat = mats[idx1 | idx2]
        if isinstance(vmin, list):
            vmin_ = vmin[c]
        else:
            vmin_ = vmin
        if isinstance(vmax, list):
            vmax_ = vmax[c]
        else:   
            vmax_ = vmax
        mat, colorbar, results = plot_pileup_mat(mat, axs[c], vmax=vmax_, vmin=vmin_, res=res, s1=s1, s2=s2,
                                                 method=center_method, center=center, show_filts=show_filts,
                                                 delta_co=delta_co)    

        all_results[u1] = results

    for c, u in enumerate(np.unique(clusters)):
        plt.sca(axs[c])
        plt.title(f'{columns_to_names.get(u, str(u)).replace(newline, " ")}', color = row_colors_dict[columns_to_names[u]])

    cbar_ax = axs[-1].inset_axes(transform = axs[c].transAxes, bounds = (1.1, 0, .05, 1))
    plt.colorbar(colorbar, cax=cbar_ax, shrink=.5)
    cbar_ax.grid(False)
    fig.add_axes(cbar_ax)
    plt.sca(cbar_ax)
    plt.yticks(fontsize=4)

    mb = mat.shape[1]*res // 2 / 1e6
    for c, a in enumerate(axs):
        cutoff = np.round(mb, 1)
        a.grid(False)
        a.set_xticks([-mb, 0, mb])
        if c%n == 0:
            a.set_yticklabels([-cutoff, "Anchor1", cutoff])
        else:
            a.set_yticklabels([])
        if c >= n*n-n:
            a.set_xticklabels([-cutoff, "Anchor2", cutoff])
            a.tick_params(labeltop = False, top = False, labelbottom = True, bottom = True)
        else:
            a.tick_params(labeltop = False, top = False, labelbottom = False, bottom = True)

        a.set_yticks([-mb, 0, mb])
        a.get_yticklabels()[1].set_rotation(90)
        a.get_yticklabels()[1].set_va('center')

    figs.append(fig)
    return figs, all_results


newline = '\n'
from compute_differential_hic_hub import make_outside_filt, get_inside_outside, test_inside_outside_baseline
def plot_pileup_mat(mats, ax, vmin=-.06, vmax=.2, cmap='gist_heat_r', res = 50_000, s1=0, s2=None, center = 2, method=False,
                    delta_co = 0, show_filts = False, log=False, cliplo=-1, cliphigh=10):
    results_dict = {}
    plt.sca(ax)
    mats = mats[:, s1:s2, s1:s2]
    mat = np.nanmean(mats, axis=0)
    if log:
        mat = np.log2(mat)
    mb = len(mat)*res / 2 / 1e6
    n = len(mat)
    cbar = ax.matshow(mat, cmap=cmap, vmin = vmin, vmax = vmax, extent = [-mb, mb, mb, -mb], )
    ax.text(.03, .05, f'n={len(mats)}', transform=ax.transAxes, fontsize=8)
    # name1 = columns_to_names[u1].replace(newline, "").replace("SE", "").replace("(", "").replace(")", "")
    # name2 = columns_to_names[u2].replace(newline, "").replace("SE", "").replace("(", "").replace(")", "")
    # plt.title(f'{name1} to {name2}')

    
    inside_filt, outside_filt = make_outside_filt(n//2, center, method=method)
    (v_middle, v_outside) = get_inside_outside(mats, inside_filt, outside_filt)
    stat, p, delta, _, _ = test_inside_outside_baseline(v_middle, v_outside, cliplo=cliplo, cliphigh=cliphigh)

    if show_filts:
        ax.matshow(inside_filt, cmap='bwr', alpha=.15, extent = [-mb, mb, mb, -mb], vmin=-1, vmax=1)
        ax.matshow(outside_filt, cmap='Purples', alpha=.1, extent = [-mb, mb, mb, -mb])


    delta = delta.round(3)
    # m = n//2
    # sl = slice(m-3, m+5+1)
    # cross = np.zeros_like(mat)
    # cross[sl, :] = 1
    # cross[:, sl] = 1
    # cross[sl, sl] = 0
    
    # X, Y = np.where(cross>0)
    # p = nonan_test(np.nanmean(mats[:, sl, sl], axis=(1, 2)),
    #                 np.nanmean(mats[:, X, Y], axis=(1))
    #                 )[1]    
    # enrich = (np.nanmean(mats[:, sl, sl], axis=(0, 1, 2)) - np.nanmean(mats[:, X, Y], axis=(0, 1))).round(3)
    pstr = format_pvalue(p)
    if np.abs(delta) < delta_co:
        pstr = 'NS'
    ax.text(.03, .95, f'MS: {delta}\np={pstr}', transform=ax.transAxes, fontsize=5, va='top')    
    
    results_dict['v_middle']  = v_middle
    results_dict['v_outside']  = v_outside
    results_dict['outside_filt'] = outside_filt
    results_dict['inside_filt'] = inside_filt

    results_dict['mat'] = mat
    # results_dict['sl'] = sl
    # results_dict['X'] = X
    # results_dict['Y'] = Y
    return mat, cbar, results_dict


# newline = '\n'
# def plot_intra_pileup_mat(mats, ax, vmin=-.06, vmax=.2, cmap='gist_heat_r', res = 50_000):
#     plt.sca(ax)
#     mat = np.nanmean(mats, axis=0)
#     mb = len(mat)*res / 2 / 1e3
#     n = len(mat)
#     cbar = ax.matshow(mat, cmap=cmap, vmin = vmin, vmax = vmax, extent = [-mb, mb, mb, -mb], )
#     ax.text(.03, .05, f'n={len(mats)}', transform=ax.transAxes, fontsize=5)
#     # name1 = columns_to_names[u1].replace(newline, "").replace("SE", "").replace("(", "").replace(")", "")
#     # name2 = columns_to_names[u2].replace(newline, "").replace("SE", "").replace("(", "").replace(")", "")
#     # plt.title(f'{name1} to {name2}')
#     m = n//2
#     sl = slice(m-3, m+5+1)
#     cross = np.zeros_like(mat)
#     cross[sl, :] = 1
#     cross[:, sl] = 1
#     cross[sl, sl] = 0
    
#     X, Y = np.where(cross>0)
#     p = nonan_test(np.nanmean(mats[:, sl, sl], axis=(1, 2)),
#                     np.nanmean(mats[:, X, Y], axis=(1))
#                     )[1]    
#     enrich = (np.nanmean(mats[:, sl, sl], axis=(0, 1, 2)) - np.nanmean(mats[:, X, Y], axis=(0, 1))).round(3)
#     ax.text(.03, .95, f'∆OE: {enrich}\np={format_pvalue(p)}', transform=ax.transAxes, fontsize=5, va='top')    
#     return cbar


from compute_differential_hic_hub import get_inside_outside
def plot_hub_violin(ind, clustdict, sep_oe_mat_treg, sep_oe_mat_tcon, outside_filt, cliplo, cliphigh, name='', 
                    fgsz=(15*mm, 15*mm), box_ylim=[-2, 2],
                    columns_to_names = {}):
    
    n = len(clustdict)
    dfs_for_boxplot = []
    pvals = []
    x_order = []
    fig, axs = init_subplots_exact(1, 1, fgsz=fgsz, dpi=300, as_list=False, xspace=1.3)
    for c, (u, goodinds) in enumerate(clustdict.items()):
        _, _, _, v_middle, v_outside = get_inside_outside(ind, goodinds, sep_oe_mat_treg, sep_oe_mat_tcon, 
                                                          outside_filt, cliplo=cliplo, cliphigh=cliphigh)
        
        df = pd.DataFrame([v_middle, v_outside], index = ['Middle', 'Outside']).T
        df = df.melt(var_name='variable', value_name='value')
        p = scipy.stats.ranksums(v_middle, v_outside, nan_policy='omit')[1]
        pvals.append(p)
        x_order.append(columns_to_names[u])
        df['x'] = columns_to_names[u]
        dfs_for_boxplot.append(df)

    plt.sca(axs)
    all_dfs = pd.concat(dfs_for_boxplot, axis=0)
    sns.boxplot(data=all_dfs, x='x', hue='variable', y='value',
                palette=['tab:purple', 'lightgray'], order=x_order,
                showfliers=False, linewidth=.2)
    for c, p in enumerate(pvals):
        if p < 1e-4:
            plt.text(c, box_ylim[1] * .9, '**', ha='center', va='center', fontsize=6,  color='black')
        elif p < .05:
            plt.text(c, box_ylim[1] * .9, '*', ha='center', va='center', fontsize=6,  color='black')

    
    ### Add asterisk to `x` which are significantly different from the outside
    plt.legend(title='', frameon=False, loc='upper left', bbox_to_anchor=(1, 1), fontsize=4)
    plt.ylabel('O/E',fontsize=4)
    plt.yticks(fontsize=4)
    plt.ylim(box_ylim)
    # plt.yticks([-1, 0, 1])
    plt.xticks(fontsize=4, rotation=20)
    add_yaxis_labels('Tcon', 'Treg', axs, fontsize=4)
    plt.title("O/E", fontsize=6)
    axs.set_axisbelow(True)
    plt.xlabel("")
    return fig






# newline = '\n'
# def plot_hic_pileup_mats(mats, idx, ax, vmin=-.06, vmax=.2, columns_to_names={}, cmap='gist_heat_r', res = 50_000):
#     plt.sca(ax)
#     mat = np.nanmean(mats[idx], axis=0)
#     mb = len(mat)*res / 2 / 1e3
#     n = len(mat)    
#     ax.matshow(mat, cmap=cmap, vmin = vmin, vmax = vmax, extent = [-mb, mb, mb, -mb])
#     ax.text(.03, .05, f'n={idx.sum()}', transform=ax.transAxes, fontsize=8)
#     m = n//2
#     sl = slice(m-3, m+5+1)
#     cross = np.zeros_like(mat)
#     cross[sl, :] = 1
#     cross[:, sl] = 1
#     cross[sl, sl] = 0
    
#     X, Y = np.where(cross>0)
#     p = nonan_test(np.nanmean(mats[idx][:, sl, sl], axis=(1, 2)),
#                     np.nanmean(mats[idx][:, X, Y], axis=(1))
#                     )[1]    
#     enrich = (np.nanmean(mats[idx][:, sl, sl], axis=(0, 1, 2)) - np.nanmean(mats[idx][:, X, Y], axis=(0, 1))).round(3)
#     ax.text(.03, .95, f'∆OE: {enrich}\np={format_pvalue(p)}', transform=ax.transAxes, fontsize=5, va='top')    
