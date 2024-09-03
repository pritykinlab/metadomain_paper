from plotting_functions import init_subplots
import matplotlib.pyplot as plt
from aux_functions import *
import numpy as np
import pandas as pd
import seaborn as sns

def foxp3_at_nearby_atac_for_joris_wei(cols, foxp3_metadata_at_all_atac_sites, foxp3_values_at_all_atac_sites, 
                                    my_tss_df,  pval_df_all, lfc_df_all, atac_bedtool, col_to_name, slopsize = 10_000,
                                    xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding'):
    tmp = foxp3_metadata_at_all_atac_sites['treg_yuri_foxp3_bw'].set_index('name')
    treg_col, tcon_col = ['treg_yuri_foxp3_bw', 'tcon_yuri_foxp3_bw']
    tss_df_set = set(my_tss_df['gene_name'])
    fig, axs = init_subplots(6, 2, space=.4, fgsz=(6, 4))
    for c, (wei_col, jor_col) in enumerate(zip(*cols)):
        wei_up = pval_df_all.index[(pval_df_all[wei_col] < .05) & (lfc_df_all[wei_col] > 0)]
        jor_up = pval_df_all.index[(pval_df_all[jor_col] < .05) & (lfc_df_all[jor_col] > 0)]
        jor_down = pval_df_all.index[(pval_df_all[jor_col] < .05) & (lfc_df_all[jor_col] < 0)]
        wei_down = pval_df_all.index[(pval_df_all[wei_col] < .05) & (lfc_df_all[wei_col] < 0)]
        
        jor_up = [x for x in jor_up if x in tss_df_set]
        wei_up = [x for x in wei_up if x in tss_df_set]
        jor_down = [x for x in jor_down if x in tss_df_set]
        wei_down = [x for x in wei_down if x in tss_df_set]  
        
        jor_up_inds = get_col(atac_bedtool.intersect(pbt.BedTool.from_dataframe(my_tss_df.set_index('gene_name').loc[wei_up]
                                            ).slop(b=slopsize, g='./annotations/chr_chromsizes',
                                                    ), u=True), -1).astype(int)
        wei_up_inds = get_col(atac_bedtool.intersect(pbt.BedTool.from_dataframe(my_tss_df.set_index('gene_name').loc[jor_up]
                                            ).slop(b=slopsize, g='./annotations/chr_chromsizes',
                                                    ), u=True), -1).astype(int)
        jor_down_inds = get_col(atac_bedtool.intersect(pbt.BedTool.from_dataframe(my_tss_df.set_index('gene_name').loc[jor_down]
                                            ).slop(b=slopsize, g='./annotations/chr_chromsizes',
                                                    ), u = True), -1).astype(int)
        wei_down_inds = get_col(atac_bedtool.intersect(pbt.BedTool.from_dataframe(my_tss_df.set_index('gene_name').loc[wei_down]
                                            ).slop(b=slopsize, g='./annotations/chr_chromsizes',
                                                    ), u = True), -1).astype(int)
        

        jor_up_vals = foxp3_values_at_all_atac_sites[treg_col][jor_up_inds] - foxp3_values_at_all_atac_sites[tcon_col][jor_up_inds]
        wei_up_vals = foxp3_values_at_all_atac_sites[treg_col][wei_up_inds] - foxp3_values_at_all_atac_sites[tcon_col][wei_up_inds]
        jor_down_vals = foxp3_values_at_all_atac_sites[treg_col][jor_down_inds] - foxp3_values_at_all_atac_sites[tcon_col][jor_down_inds]
        wei_down_vals = foxp3_values_at_all_atac_sites[treg_col][wei_down_inds] - foxp3_values_at_all_atac_sites[tcon_col][wei_down_inds]
        baseline = foxp3_values_at_all_atac_sites[treg_col] - foxp3_values_at_all_atac_sites[tcon_col]
        xs = np.linspace(-500, 500, jor_up_vals.shape[1])
        plt.sca(axs[c])
        plt.plot(xs, np.nanmean(jor_up_vals, axis=0), color='#FF9999', label=f'KO \n(n={len(jor_up_vals)})')
        plt.plot(xs, np.nanmean(wei_up_vals, axis=0), color='#FF6633', label=f'Tir1 \n(n={len(wei_up_vals)})')
        plt.plot(xs, np.nanmean(baseline, axis=0), color='black')
        plt.legend(loc='upper left')
        
    #     p_ks = scipy.stats.ks_2samp(jor_up_vals[:, 500], wei_up_vals[:, 500])[1]
    #     p_t = scipy.stats.ttest_ind(jor_up_vals[:, 500], wei_up_vals[:, 500])[1]
        p_t_baseline = scipy.stats.ttest_ind(baseline[:, 500], wei_up_vals[:, 500])[1]
        p_t_baseline2 = scipy.stats.ttest_ind(baseline[:, 500], jor_up_vals[:, 500])[1]
    #     plt.text(.98, 1-.02, f'Jor/Wei_ks={format_pvalue(p_ks)}', ha='right', va='top', transform=plt.gca().transAxes)
    #     plt.text(.98, 1-.12, f'Jor/Wei_t={format_pvalue(p_t)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.02, f'Wei/Base={format_pvalue(p_t_baseline)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.12, f'Jor/Base={format_pvalue(p_t_baseline2)}', ha='right', va='top', transform=plt.gca().transAxes)

        
        plt.title(col_to_name[wei_col] + " Up")
        
        plt.sca(axs[c+3])
        plt.plot(xs, np.nanmean(jor_down_vals, axis=0), color='#9999FF', label=f'KO \n(n={len(jor_down_vals)})')
        plt.plot(xs, np.nanmean(wei_down_vals, axis=0), color='#66CCFF', label=f'Tir1 \n(n={len(wei_down_vals)})')
        plt.plot(xs, np.nanmean(baseline, axis=0), color='black')
        plt.title(col_to_name[wei_col] + " Down")
        plt.legend(loc='upper left')
        
        p_t_baseline = scipy.stats.ttest_ind(baseline[:, 500], wei_down_vals[:, 500])[1]
        p_t_baseline2 = scipy.stats.ttest_ind(baseline[:, 500], jor_down_vals[:, 500])[1]
    #     plt.text(.98, 1-.02, f'Jor/Wei_ks={format_pvalue(p_ks)}', ha='right', va='top', transform=plt.gca().transAxes)
    #     plt.text(.98, 1-.12, f'Jor/Wei_t={format_pvalue(p_t)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.02, f'Wei/Base={format_pvalue(p_t_baseline)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.12, f'Jor/Base={format_pvalue(p_t_baseline2)}', ha='right', va='top', transform=plt.gca().transAxes)
        for ax in axs:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            


def foxp3_at_nearest_atac_peaks(cols, foxp3_metadata_at_all_atac_sites, foxp3_values_at_all_atac_sites, 
                                    pval_df_all, lfc_df_all, atac_bedtool, col_to_name, my_atac_meta_new, 
                                    foxp3_count_set, 
                                    xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                    pco = .05):
    tmp = foxp3_metadata_at_all_atac_sites['treg_yuri_foxp3_bw']#.set_index('name')
    treg_col, tcon_col = ['treg_yuri_foxp3_bw', 'tcon_yuri_foxp3_bw']
    fig, axs = init_subplots(6, 2, space=.4, fgsz=(6, 4))
    for c, (wei_col, jor_col) in enumerate(zip(*cols)):
        wei_up = pval_df_all.index[(pval_df_all[wei_col] < pco) & (lfc_df_all[wei_col] > 0)]
        jor_up = pval_df_all.index[(pval_df_all[jor_col] < pco) & (lfc_df_all[jor_col] > 0)]
        jor_down = pval_df_all.index[(pval_df_all[jor_col] < pco) & (lfc_df_all[jor_col] < 0)]
        wei_down = pval_df_all.index[(pval_df_all[wei_col] < pco) & (lfc_df_all[wei_col] < 0)]
        
        jor_up = [x for x in jor_up if x in foxp3_count_set]
        wei_up = [x for x in wei_up if x in foxp3_count_set]
        jor_down = [x for x in jor_down if x in foxp3_count_set]
        wei_down = [x for x in wei_down if x in foxp3_count_set]  
        jor_up_inds = tmp.index.isin(my_atac_meta_new['index'][my_atac_meta_new[0].isin(jor_up)])
        wei_up_inds = tmp.index.isin(my_atac_meta_new['index'][my_atac_meta_new[0].isin(wei_up)])
        jor_down_inds = tmp.index.isin(my_atac_meta_new['index'][my_atac_meta_new[0].isin(jor_down)])
        wei_down_inds = tmp.index.isin(my_atac_meta_new['index'][my_atac_meta_new[0].isin(wei_down)])

        jor_up_vals = foxp3_values_at_all_atac_sites[treg_col][jor_up_inds] - foxp3_values_at_all_atac_sites[tcon_col][jor_up_inds]

        wei_up_vals = foxp3_values_at_all_atac_sites[treg_col][wei_up_inds] - foxp3_values_at_all_atac_sites[tcon_col][wei_up_inds]
        jor_down_vals = foxp3_values_at_all_atac_sites[treg_col][jor_down_inds] - foxp3_values_at_all_atac_sites[tcon_col][jor_down_inds]
        wei_down_vals = foxp3_values_at_all_atac_sites[treg_col][wei_down_inds] - foxp3_values_at_all_atac_sites[tcon_col][wei_down_inds]
        baseline = foxp3_values_at_all_atac_sites[treg_col] - foxp3_values_at_all_atac_sites[tcon_col]
        xs = np.linspace(-500, 500, jor_up_vals.shape[1])
        plt.sca(axs[c])
        plt.plot(xs, np.nanmean(jor_up_vals, axis=0), color='#FF9999', label=f'KO \n(n={len(jor_up_vals)})')
        plt.plot(xs, np.nanmean(wei_up_vals, axis=0), color='#FF6633', label=f'Tir1 \n(n={len(wei_up_vals)})')
        plt.plot(xs, np.nanmean(baseline, axis=0), color='black')
        plt.legend(loc='upper left')
        
    #     p_ks = scipy.stats.ks_2samp(jor_up_vals[:, 500], wei_up_vals[:, 500])[1]
    #     p_t = scipy.stats.ttest_ind(jor_up_vals[:, 500], wei_up_vals[:, 500])[1]
        p_t_baseline = scipy.stats.ttest_ind(baseline[:, 500], wei_up_vals[:, 500])[1]
        p_t_baseline2 = scipy.stats.ttest_ind(baseline[:, 500], jor_up_vals[:, 500])[1]
    #     plt.text(.98, 1-.02, f'Jor/Wei_ks={format_pvalue(p_ks)}', ha='right', va='top', transform=plt.gca().transAxes)
    #     plt.text(.98, 1-.12, f'Jor/Wei_t={format_pvalue(p_t)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.02, f'Wei/Base={format_pvalue(p_t_baseline)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.12, f'Jor/Base={format_pvalue(p_t_baseline2)}', ha='right', va='top', transform=plt.gca().transAxes)

        
        plt.title(col_to_name[wei_col] + " Up")
        
        plt.sca(axs[c+3])
        plt.plot(xs, np.nanmean(jor_down_vals, axis=0), color='#9999FF', label=f'KO \n(n={len(jor_down_vals)})')
        plt.plot(xs, np.nanmean(wei_down_vals, axis=0), color='#66CCFF', label=f'Tir1 \n(n={len(wei_down_vals)})')
        plt.plot(xs, np.nanmean(baseline, axis=0), color='black')
        plt.title(col_to_name[wei_col] + " Down")
        plt.legend(loc='upper left')
        
        p_t_baseline = scipy.stats.ttest_ind(baseline[:, 500], wei_down_vals[:, 500])[1]
        p_t_baseline2 = scipy.stats.ttest_ind(baseline[:, 500], jor_down_vals[:, 500])[1]
    #     plt.text(.98, 1-.02, f'Jor/Wei_ks={format_pvalue(p_ks)}', ha='right', va='top', transform=plt.gca().transAxes)
    #     plt.text(.98, 1-.12, f'Jor/Wei_t={format_pvalue(p_t)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.02, f'Wei/Base={format_pvalue(p_t_baseline)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.12, f'Jor/Base={format_pvalue(p_t_baseline2)}', ha='right', va='top', transform=plt.gca().transAxes)
        for ax in axs:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)


def n_tfs_at_degs(cols,  pval_df_all, lfc_df_all, basemean_df_all,
                                    col_to_name, tf_df,
                                    xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding', pco=0.05,
                                    basemean_co = 20):
    n = len(cols[0])
    fig, axs = init_subplots(n, 1, space=.6, fgsz=(8, 4))

    dict_of_vals = {}
    for c, (wei_col, jor_col) in enumerate(zip(*cols)):
        plt.sca(axs[c])
        basemean_idx = (basemean_df_all[wei_col] > basemean_co) & (basemean_df_all[jor_col] > basemean_co)
        high_basemean_genes = basemean_df_all.index[basemean_idx]
        pval_df, lfc_df = pval_df_all.loc[high_basemean_genes], lfc_df_all.loc[high_basemean_genes]

        wei_up = pval_df.index[(pval_df[wei_col] < pco) & (lfc_df[wei_col] > 0)]
        jor_up = pval_df.index[(pval_df[jor_col] < pco) & (lfc_df[jor_col] > 0)]
        both_up = wei_up.intersection(jor_up)

        jor_down = pval_df.index[(pval_df[jor_col] < pco) & (lfc_df[jor_col] < 0)]
        wei_down = pval_df.index[(pval_df[wei_col] < pco) & (lfc_df[wei_col] < 0)]
        both_down = wei_down.intersection(jor_down)
        baseline_genes = pval_df.index

        jor_up_frac = np.isin(jor_up, tf_df.index).mean()
        wei_up_frac = np.isin(wei_up, tf_df.index).mean()
        both_up_frac = np.isin(both_up, tf_df.index).mean()
        jor_down_frac = np.isin(jor_down, tf_df.index).mean()
        wei_down_frac = np.isin(wei_down, tf_df.index).mean()
        both_down_frac = np.isin(both_down, tf_df.index).mean()
        baseline_frac = np.isin(baseline_genes, tf_df.index).mean()

        jor_up_tf_count = np.isin(jor_up, tf_df.index).sum()
        jor_up_no_tf_count = (~np.isin(jor_up, tf_df.index)).sum()
        wei_up_tf_count = np.isin(wei_up, tf_df.index).sum()
        wei_up_no_tf_count = (~np.isin(wei_up, tf_df.index)).sum()
        both_up_tf_count = np.isin(both_up, tf_df.index).sum()
        both_up_no_tf_count = (~np.isin(both_up, tf_df.index)).sum()
        jor_down_tf_count = np.isin(jor_down, tf_df.index).sum()
        jor_down_no_tf_count = (~np.isin(jor_down, tf_df.index)).sum()
        wei_down_tf_count = np.isin(wei_down, tf_df.index).sum()
        wei_down_no_tf_count = (~np.isin(wei_down, tf_df.index)).sum()
        both_down_tf_count = np.isin(both_down, tf_df.index).sum()
        both_down_no_tf_count = (~np.isin(both_down, tf_df.index)).sum()
        baseline_tf_count = np.isin(baseline_genes, tf_df.index).sum()
        baseline_no_tf_count = (~np.isin(baseline_genes, tf_df.index)).sum()

        _, pval_jor_up = scipy.stats.fisher_exact([[jor_up_tf_count, jor_up_no_tf_count], [baseline_tf_count, baseline_no_tf_count]])
        _, pval_wei_up = scipy.stats.fisher_exact([[wei_up_tf_count, wei_up_no_tf_count], [baseline_tf_count, baseline_no_tf_count]])
        _, pval_both_up = scipy.stats.fisher_exact([[both_up_tf_count, both_up_no_tf_count], [baseline_tf_count, baseline_no_tf_count]])
        _, pval_jor_down = scipy.stats.fisher_exact([[jor_down_tf_count, jor_down_no_tf_count], [baseline_tf_count, baseline_no_tf_count]])
        _, pval_wei_down = scipy.stats.fisher_exact([[wei_down_tf_count, wei_down_no_tf_count], [baseline_tf_count, baseline_no_tf_count]])
        _, pval_both_down = scipy.stats.fisher_exact([[both_down_tf_count, both_down_no_tf_count], [baseline_tf_count, baseline_no_tf_count]])
        
        plt.bar('baseline', baseline_frac, color='lightgray', zorder=3)
        plt.bar('jor_up', jor_up_frac, color='#FF9999', zorder=3, label=f'jor_up: {format_pvalue(pval_jor_up)}')
        plt.bar('wei_up', wei_up_frac, color='#FF6633', zorder=3, label=f'wei_up: {format_pvalue(pval_wei_up)}')
        plt.bar('both_up', both_up_frac, color='#FF0000', zorder=3, label=f'both_up: {format_pvalue(pval_both_up)}')
        plt.bar('jor_down', jor_down_frac, color='#9999FF', zorder=3, label=f'jor_down: {format_pvalue(pval_jor_down)}')
        plt.bar('wei_down', wei_down_frac, color='#66CCFF', zorder=3, label=f'wei_down: {format_pvalue(pval_wei_down)}')
        plt.bar('both_down', both_down_frac, color='#0000FF', zorder=3, label=f'both_down: {format_pvalue(pval_both_down)}')

        
        plt.ylabel('Fraction DEGs which are TFs')
        plt.title(col_to_name[wei_col].split(".")[0])
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)



def n_tfs_at_degs_over_time(cols,  pval_df_all, lfc_df_all, basemean_df_all,
                                    col_to_name, tf_df,
                                    xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding', pco=0.05,
                                    basemean_co = 20):
    fig, axs = init_subplots(1, 1, space=.6, fgsz=(8, 4))

    dict_of_vals = {}
    n = len(cols)
    colors = sns.color_palette("coolwarm", n_colors=n*2)

    for c, (col) in enumerate(cols):
        basemean_idx = (basemean_df_all[col] > basemean_co) & (basemean_df_all[col] > basemean_co)
        high_basemean_genes = basemean_df_all.index[basemean_idx]
        pval_df, lfc_df = pval_df_all.loc[high_basemean_genes], lfc_df_all.loc[high_basemean_genes]

        up = pval_df.index[(pval_df[col] < pco) & (lfc_df[col] > 0)]
        down = pval_df.index[(pval_df[col] < pco) & (lfc_df[col] < 0)]
        baseline_genes = pval_df.index
        up_frac = np.isin(up, tf_df.index).mean()
        down_frac = np.isin(down, tf_df.index).mean()
        baseline_frac = np.isin(baseline_genes, tf_df.index).mean()

        plt.bar(col_to_name[col].split(" ")[1] + ' up', up_frac, color=colors[c], zorder=3)
        plt.bar(col_to_name[col].split(" ")[1] + ' down', down_frac, color=colors[2*n-1-c], zorder=3)

    plt.bar('baseline', baseline_frac, color='lightgray', zorder=3)
    plt.ylabel('Fraction DEGs which are TFs')



def n_foxp3_peaks_near_degs(cols,  pval_df_all, lfc_df_all, basemean_df_all,
                                    col_to_name, foxp3_count_set, full_foxp3_counts,
                                    xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding', pco=0.05,
                                    basemean_co = 20,
                                    label1='KO', label2='Tir1', label3='Both'):
    n = len(cols[0])
    fig, axs = init_subplots(n*2, 2, space=.4, fgsz=(7, 4))

    dict_of_vals = {}
    for c, (wei_col, jor_col) in enumerate(zip(*cols)):
        basemean_idx = (basemean_df_all[wei_col] > basemean_co) & (basemean_df_all[jor_col] > basemean_co)
        high_basemean_genes = basemean_df_all.index[basemean_idx]
        pval_df, lfc_df = pval_df_all.loc[high_basemean_genes], lfc_df_all.loc[high_basemean_genes]

        wei_up = pval_df.index[(pval_df[wei_col] < pco) & (lfc_df[wei_col] > 0)]
        jor_up = pval_df.index[(pval_df[jor_col] < pco) & (lfc_df[jor_col] > 0)]
        both_up = wei_up.intersection(jor_up)

        jor_down = pval_df.index[(pval_df[jor_col] < pco) & (lfc_df[jor_col] < 0)]
        wei_down = pval_df.index[(pval_df[wei_col] < pco) & (lfc_df[wei_col] < 0)]
        both_down = wei_down.intersection(jor_down)

        wei_up_counts_to_save = [(x, full_foxp3_counts.loc[x]) if x in foxp3_count_set else (x, "LOW_EXPR") for x in wei_up]
        wei_down_counts_to_save = [(x, full_foxp3_counts.loc[x]) if x in foxp3_count_set else (x, "LOW_EXPR") for x in wei_down]

        jor_up = [x for x in jor_up if x in foxp3_count_set]
        wei_up = [x for x in wei_up if x in foxp3_count_set]
        both_up = [x for x in both_up if x in foxp3_count_set]
        jor_down = [x for x in jor_down if x in foxp3_count_set]
        wei_down = [x for x in wei_down if x in foxp3_count_set]  
        both_down = [x for x in both_down if x in foxp3_count_set]
        
        dict_of_vals[col_to_name[wei_col] + " Up"] = wei_up_counts_to_save
        dict_of_vals[col_to_name[wei_col] + " Down"] = wei_down_counts_to_save

        jor_up_foxp3_counts = full_foxp3_counts.loc[jor_up]
        wei_up_foxp3_counts = full_foxp3_counts.loc[wei_up]
        both_up_foxp3_counts = full_foxp3_counts.loc[both_up]

        jor_down_foxp3_counts = full_foxp3_counts.loc[jor_down]
        wei_down_foxp3_counts = full_foxp3_counts.loc[wei_down]
        both_down_foxp3_counts = full_foxp3_counts.loc[both_down]

        baseline = full_foxp3_counts.loc[[x for x in pval_df_all.index if x in foxp3_count_set]]
        
        plt.sca(axs[c])



        v = plot_foxp3_peaks_near_degs_as_bar(jor_up, col_to_name, label1, foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = '#FF9999',
                                           linewidth = 3
                               )
        v = plot_foxp3_peaks_near_degs_as_bar(wei_up, col_to_name, label2, foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = '#FF6633',
                                           linewidth = 3
                               )
        v = plot_foxp3_peaks_near_degs_as_bar(both_up, col_to_name, label3, foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = '#FF0000',
                                           linewidth = 3
                               )
        v = plot_foxp3_peaks_near_degs_as_bar(foxp3_count_set, col_to_name, 'degs', foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = 'lightgray',
                                           linewidth = 3, plot_baseline=True
                               )

        # sns.ecdfplot(jor_up_foxp3_counts, color='#FF9999', label=f'KO Up (n={len(jor_up_foxp3_counts)})')
        # # sns.ecdfplot(wei_up_foxp3_counts, color='#FF6633', label=f'Tir1 Up (n={len(wei_up_foxp3_counts)})')
        # sns.ecdfplot(both_up_foxp3_counts, color='#FF0000', label=f'Both Up (n={len(both_up_foxp3_counts)})')
        # sns.ecdfplot(baseline, color='black', label=f'Baseline')
        # if len(wei_up_foxp3_counts) == 0:
            # continue
        
        p_wei = scipy.stats.ranksums(wei_up_foxp3_counts, baseline)[1]
        p_jor = scipy.stats.ranksums(jor_up_foxp3_counts, baseline)[1]
        p_both = scipy.stats.ranksums(both_up_foxp3_counts, baseline)[1]

        # plt.text(.98, .02, f'wei={format_pvalue(p_wei)}', ha='right', va='bottom', transform=plt.gca().transAxes)
        # plt.text(.98, .12, f'jor={format_pvalue(p_jor)}', ha='right', va='bottom', transform=plt.gca().transAxes)
        # plt.text(.98, .22, f'both={format_pvalue(p_both)}', ha='right', va='bottom', transform=plt.gca().transAxes)

        plt.title(f"{col_to_name[wei_col]} Up")
        plt.legend(loc = 'upper left', 
                    # bbox_to_anchor=(1., 0), 
                    frameon=False)
        
        plt.sca(axs[c+n])

        v = plot_foxp3_peaks_near_degs_as_bar(jor_down, col_to_name, label1, foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = '#9999FF',
                                           linewidth = 3
                               )
        v = plot_foxp3_peaks_near_degs_as_bar(wei_down, col_to_name, label2, foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = '#66CCFF',
                                           linewidth = 3
                               )
        v = plot_foxp3_peaks_near_degs_as_bar(both_down, col_to_name, label3, foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = '#0000FF',
                                           linewidth = 3
                               )
        v = plot_foxp3_peaks_near_degs_as_bar(foxp3_count_set, col_to_name, 'degs', foxp3_count_set, full_foxp3_counts,
                                            xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                            color = 'lightgray',
                                           linewidth = 3, plot_baseline=True
                               )

        # sns.ecdfplot(jor_down_foxp3_counts, color='#9999FF', label=f'KO Down (n={len(jor_down_foxp3_counts)})')
        # sns.ecdfplot(wei_down_foxp3_counts, color='#66CCFF', label=f'Tir1 Down (n={len(wei_down_foxp3_counts)})')
        # sns.ecdfplot(both_down_foxp3_counts, color='#0000FF', label=f'Both Down (n={len(both_down_foxp3_counts)})')
        # sns.ecdfplot(baseline, color='black', label=f'Baseline')
        plt.title(f"{col_to_name[wei_col]} Down")
        p_wei = scipy.stats.ranksums(wei_down_foxp3_counts, baseline)[1]
        p_jor = scipy.stats.ranksums(jor_down_foxp3_counts, baseline)[1]
        p_both = scipy.stats.ranksums(both_down_foxp3_counts, baseline)[1]

        # plt.text(.98, .02, f'wei={format_pvalue(p_wei)}', ha='right', va='bottom', transform=plt.gca().transAxes)
        # plt.text(.98, .12, f'jor={format_pvalue(p_jor)}', ha='right', va='bottom', transform=plt.gca().transAxes)
        # plt.text(.98, .22, f'both={format_pvalue(p_both)}', ha='right', va='bottom', transform=plt.gca().transAxes)
        plt.legend(loc = 'lower left', bbox_to_anchor=(1, 0), frameon=False)
        for ax in axs:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            # ax.set_xlim([-.2, 10])
    for ax in axs:
        plt.sca(ax)
        plt.xticks([-1, 0, 1, 2])
        z = list(ax.get_xticklabels())
        z[0].set_text("All")
        ax.set_xticklabels(z)
    return fig, dict_of_vals


def plot_foxp3_peaks_near_degs(degs, col_to_name, foxp3_count_set, full_foxp3_counts, ax=None, color=None,
                                    xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding', linewidth=2):
    ax = plt.gca() if ax is None else ax
    # fig, axs = init_subplots(6, 2, space=.4, fgsz=(6, 4))
    degs = [x for x in degs if x in foxp3_count_set]
    if len(degs)==0:
        return 0
    baseline = full_foxp3_counts.loc[[x for x in foxp3_count_set]]
    wei_up_foxp3_counts = full_foxp3_counts.loc[degs]
    p_wei = scipy.stats.ranksums(wei_up_foxp3_counts, baseline)[1]
    
    # plt.sca(axs[c])
    if p_wei < .05:
        sns.ecdfplot(wei_up_foxp3_counts, label=f'(p={format_pvalue(p_wei)})', c = color, linewidth=linewidth)
    else:
        sns.ecdfplot(wei_up_foxp3_counts, c = color, linewidth=linewidth)
    # sns.ecdfplot(wei_up_foxp3_counts, color='#FF6633', label=f'Tir1 Up (n={len(wei_up_foxp3_counts)})')
    sns.ecdfplot(baseline, color='black', linewidth=linewidth)
    # if len(wei_up_foxp3_counts) == 0:
        # continue
    # p_jor = scipy.stats.ks_2samp(jor_up_foxp3_counts, baseline)[1]
    # plt.text(.98, .02, f'wei={format_pvalue(p_wei)}', ha='right', va='bottom', transform=plt.gca().transAxes)
    # plt.text(.98, .12, f'jor={format_pvalue(p_jor)}', ha='right', va='bottom', transform=plt.gca().transAxes)

    # plt.title(col_to_name[wei_col] + " Up")
    plt.legend(loc = 'lower right', bbox_to_anchor=(1.02, .1), frameon=False)
    
    # plt.sca(axs[c+3])
    # sns.ecdfplot(jor_down_foxp3_counts, color='#9999FF', label=f'KO Down (n={len(jor_down_foxp3_counts)})')
    # sns.ecdfplot(wei_down_foxp3_counts, color='#66CCFF', label=f'Tir1 Down (n={len(wei_down_foxp3_counts)})')
    # sns.ecdfplot(baseline, color='black', label=f'Baseline')
    # plt.title(col_to_name[wei_col] + " Down")
    # p_wei = scipy.stats.ks_2samp(wei_down_foxp3_counts, baseline)[1]
    # p_jor = scipy.stats.ks_2samp(jor_down_foxp3_counts, baseline)[1]
    # plt.text(.98, .02, f'wei={format_pvalue(p_wei)}', ha='right', va='bottom', transform=plt.gca().transAxes)
    # plt.text(.98, .12, f'jor={format_pvalue(p_jor)}', ha='right', va='bottom', transform=plt.gca().transAxes)
    # plt.legend(loc = 'lower right', bbox_to_anchor=(1.02, .1), frameon=False)
    # for ax in axs:
    #     ax.set_xlabel(xlabel)
    #     ax.set_ylabel(ylabel)
    #     ax.set_xlim([-.2, 10])

def plot_foxp3_peaks_near_degs_as_bar(degs, col_to_name, name, foxp3_count_set, full_foxp3_counts, ax=None, color=None,
                                    xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding', linewidth=2, plot_baseline=False):
    ax = plt.gca() if ax is None else ax
    # fig, axs = init_subplots(6, 2, space=.4, fgsz=(6, 4))
    degs = [x for x in degs if x in foxp3_count_set]
    if len(degs)==0:
        return 0
    baseline = full_foxp3_counts.loc[[x for x in foxp3_count_set]]
    wei_up_foxp3_counts = full_foxp3_counts.loc[degs]
    p_wei = scipy.stats.ranksums(wei_up_foxp3_counts, baseline)[1]
    
    # plt.sca(axs[c])
    # if p_wei < .05:
    if plot_baseline==True:
        sem = baseline.sem()  # Standard deviation
        plt.bar(-1, baseline.mean(), color='lightgray', linewidth=linewidth, zorder=3)
        plt.errorbar(-1, baseline.mean(), yerr=sem, fmt='o', color='black', capsize=5, zorder=4, markersize=0)
    else:
        sem = wei_up_foxp3_counts.sem()  # Standard deviation
        plt.bar(name, wei_up_foxp3_counts.mean(), label=f'p={format_pvalue(p_wei)}; n={len(degs)}', color = color, linewidth=linewidth,
                    zorder=3)
        plt.errorbar(name, wei_up_foxp3_counts.mean(), yerr=sem, fmt='o', color='black', capsize=5, zorder=4, markersize=0)

    plt.legend(loc = 'lower left', bbox_to_anchor=(1.02, .1), frameon=False)
    return wei_up_foxp3_counts
    # plt.sca(axs[c+3])
    # sns.ecdfplot(jor_down_foxp3_counts, color='#9999FF', label=f'KO Down (n={len(jor_down_foxp3_counts)})')
    # sns.ecdfplot(wei_down_foxp3_counts, color='#66CCFF', label=f'Tir1 Down (n={len(wei_down_foxp3_counts)})')
    # sns.ecdfplot(baseline, color='black', label=f'Baseline')
    # plt.title(col_to_name[wei_col] + " Down")
    # p_wei = scipy.stats.ks_2samp(wei_down_foxp3_counts, baseline)[1]
    # p_jor = scipy.stats.ks_2samp(jor_down_foxp3_counts, baseline)[1]
    # plt.text(.98, .02, f'wei={format_pvalue(p_wei)}', ha='right', va='bottom', transform=plt.gca().transAxes)
    # plt.text(.98, .12, f'jor={format_pvalue(p_jor)}', ha='right', va='bottom', transform=plt.gca().transAxes)
    # plt.legend(loc = 'lower right', bbox_to_anchor=(1.02, .1), frameon=False)
    # for ax in axs:
    #     ax.set_xlabel(xlabel)
    #     ax.set_ylabel(ylabel)
    #     ax.set_xlim([-.2, 10])




def get_total_gene_pileup(geneset, foxp3_values_at_all_atac_sites, tmp, my_atac_meta_new, treg_col, tcon_col):
    total_gene_pileup = []
    for gene in geneset:
        gene_atac_idx = my_atac_meta_new[0] == gene
        gene_inds = my_atac_meta_new['index'][gene_atac_idx].astype(int)
        gene_idx = tmp.index.isin(gene_inds)
        gene_vals = foxp3_values_at_all_atac_sites[treg_col][gene_idx] - foxp3_values_at_all_atac_sites[tcon_col][gene_idx]
        assert len(gene_vals.shape) == 2

        one_gene_pileup = np.nanmean(gene_vals, axis=0)
        total_gene_pileup.append(one_gene_pileup)

    return np.asarray(total_gene_pileup)


def pileup_stat_test(p1, p2, sl = slice(490, 510), test=scipy.stats.ranksums):
    p1 = np.nanmean(p1[:, sl], axis=1)
    p2 = np.nanmean(p2[:, sl], axis=1)
    bad1 = np.isnan(p1)
    bad2 = np.isnan(p2)
    p1 = p1[~bad1]
    p2 = p2[~bad2]
    return test(p1, p2)

def sum_foxp3_at_nearest_atac_peaks(cols, foxp3_metadata_at_all_atac_sites, foxp3_values_at_all_atac_sites, 
                                    pval_df_all, lfc_df_all, basemean_df_all, atac_bedtool, col_to_name, my_atac_meta_new, 
                                    foxp3_count_set, xlabel = "ATAC peaks near DEGs", ylabel='Foxp3 Binding',
                                    pco = .05, basemean_co = 20):
    tmp = foxp3_metadata_at_all_atac_sites['treg_yuri_foxp3_bw']#.set_index('name')
    treg_col, tcon_col = ['treg_yuri_foxp3_bw', 'tcon_yuri_foxp3_bw']
    fig, axs = init_subplots(6, 2, space=.4, fgsz=(6, 4))
    
    baseline_geneset = sorted(my_atac_meta_new[0].unique())

    baseline = get_total_gene_pileup(baseline_geneset, foxp3_values_at_all_atac_sites, tmp, my_atac_meta_new, treg_col, tcon_col)
    # baseline = foxp3_values_at_all_atac_sites[treg_col] - foxp3_values_at_all_atac_sites[tcon_col]
    for c, (wei_col, jor_col) in enumerate(zip(*cols)):
        basemean_idx = (basemean_df_all[wei_col] > basemean_co) & (basemean_df_all[jor_col] > basemean_co)
        high_basemean_genes = basemean_df_all.index[basemean_idx]

        pval_df, lfc_df = pval_df_all.loc[high_basemean_genes], lfc_df_all.loc[high_basemean_genes]

        wei_up = pval_df.index[(pval_df[wei_col] < pco) & (lfc_df[wei_col] > 0)]
        jor_up = pval_df.index[(pval_df[jor_col] < pco) & (lfc_df[jor_col] > 0)]
        jor_down = pval_df.index[(pval_df[jor_col] < pco) & (lfc_df[jor_col] < 0)]
        wei_down = pval_df.index[(pval_df[wei_col] < pco) & (lfc_df[wei_col] < 0)]
        
        jor_up = [x for x in jor_up if x in foxp3_count_set]
        wei_up = [x for x in wei_up if x in foxp3_count_set]
        jor_down = [x for x in jor_down if x in foxp3_count_set]
        wei_down = [x for x in wei_down if x in foxp3_count_set]  

        jor_up_vals = get_total_gene_pileup(jor_up, foxp3_values_at_all_atac_sites, tmp, my_atac_meta_new, treg_col, tcon_col)
        wei_up_vals = get_total_gene_pileup(wei_up, foxp3_values_at_all_atac_sites, tmp, my_atac_meta_new, treg_col, tcon_col)
        jor_down_vals = get_total_gene_pileup(jor_down, foxp3_values_at_all_atac_sites, tmp, my_atac_meta_new, treg_col, tcon_col)
        wei_down_vals = get_total_gene_pileup(wei_down, foxp3_values_at_all_atac_sites, tmp, my_atac_meta_new, treg_col, tcon_col)

        xs = np.linspace(-500, 500, jor_up_vals.shape[1])
        plt.sca(axs[c])
        plt.plot(xs, np.nanmean(jor_up_vals, axis=0), color='#FF9999', label=f'KO \n(n={len(jor_up_vals)})')
        plt.plot(xs, np.nanmean(wei_up_vals, axis=0), color='#FF6633', label=f'Tir1 \n(n={len(wei_up_vals)})')
        plt.plot(xs, np.nanmean(baseline, axis=0), color='black')
        plt.legend(loc='upper left')
        
    #     p_ks = scipy.stats.ks_2samp(jor_up_vals[:, 500], wei_up_vals[:, 500])[1]
    #     p_t = scipy.stats.ttest_ind(jor_up_vals[:, 500], wei_up_vals[:, 500])[1]
        vals = range(450, 550)
        # print(wei_up_vals.shape)
        p_t_baseline = pileup_stat_test(baseline, wei_up_vals)[1]
        p_t_baseline2 = pileup_stat_test(baseline, jor_up_vals)[1]
    #     plt.text(.98, 1-.02, f'Jor/Wei_ks={format_pvalue(p_ks)}', ha='right', va='top', transform=plt.gca().transAxes)
    #     plt.text(.98, 1-.12, f'Jor/Wei_t={format_pvalue(p_t)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.02, f'Wei/Base={format_pvalue(p_t_baseline)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.12, f'Jor/Base={format_pvalue(p_t_baseline2)}', ha='right', va='top', transform=plt.gca().transAxes)

        
        plt.title(col_to_name[wei_col] + " Up")
        
        plt.sca(axs[c+3])
        plt.plot(xs, np.nanmean(jor_down_vals, axis=0), color='#9999FF', label=f'KO \n(n={len(jor_down_vals)})')
        plt.plot(xs, np.nanmean(wei_down_vals, axis=0), color='#66CCFF', label=f'Tir1 \n(n={len(wei_down_vals)})')
        plt.plot(xs, np.nanmean(baseline, axis=0), color='black')
        plt.title(col_to_name[wei_col] + " Down")
        plt.legend(loc='upper left')
        
        p_t_baseline = pileup_stat_test(baseline, wei_up_vals)[1]
        p_t_baseline2 = pileup_stat_test(baseline, wei_up_vals)[1]
    #     plt.text(.98, 1-.02, f'Jor/Wei_ks={format_pvalue(p_ks)}', ha='right', va='top', transform=plt.gca().transAxes)
    #     plt.text(.98, 1-.12, f'Jor/Wei_t={format_pvalue(p_t)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.02, f'Wei/Base={format_pvalue(p_t_baseline)}', ha='right', va='top', transform=plt.gca().transAxes)
        plt.text(.98, 1-.12, f'Jor/Base={format_pvalue(p_t_baseline2)}', ha='right', va='top', transform=plt.gca().transAxes)
        for ax in axs:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        return baseline, wei_down_vals