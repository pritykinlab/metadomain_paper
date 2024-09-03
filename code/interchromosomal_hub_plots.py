from megaloop_plotting_functions import *
from plotting_functions import *
import itertools

def plot_repressive_hub(all_ind_to_region, all_gene_names, get_treg_mat, get_tcon_mat, treg_5kb, tcon_5kb, deltas, h3k27me3_bbdict):
    all_gene_names_except_these = all_gene_names - set(['Prex1', 'C1qa',  'Epha8',  'Megf11', 'Rbfox3',  'Kif26a', 'Kcnj6', 'Krt84', 'Cacna1h'])
    d = 600_000
    res = 100_000
    inds = [1447, 2699, 5344, 6582, 7045, 8466, 8856, 8971]
    deltas = [0, -100_000, 100_000, -100_000, -150_000, -60_000, -100_000, -100_000]
    deltas = dict(zip(inds, deltas))

    n = len(inds)
    label = ''
    fig, axs = init_subplots_exact(n*n, n, fgsz = (2*30*mm, 2*30*mm), space=1.1)
    for c, (i1, i2) in enumerate(itertools.product(inds, inds)):
        if (i1 != i2) and (all_ind_to_region[i1][0] == all_ind_to_region[i2][0]):
            vmin = 1e-5
            vmax = 5e-5
        elif (i1 != i2):
            vmin = 1e-5
            vmax = 6e-4
            print(i1, i2)
        else:
            vmin = 1e-4
            vmax = 5e-1
        dL = deltas[i1]
        dR = deltas[i2]
        d_LR, d_LL, d_RR, d_RL = d-dL, d+dL, d-dR, d+dR
        if i1 == i2:
            ax = axs[c]
            _, ((ax1, ax2), bw_ax, place1, place2), ax = make_scoping_plot_for_fig1_near_diag(treg_5kb, tcon_5kb, i1, i2, get_treg_mat,
                                                                                        get_tcon_mat,
                                                                                    label, resolution=10_000, 
                                                                                    d_LR=d_LR,
                                                                                    d_LL=d_LL,
                                                                                    d_RR=d_RR,
                                                                                    d_RL=d_RL,
                                        ignore_set = all_gene_names_except_these, 
                                        vmin=vmin, vmax=vmax, 
                                        vert_bbdict={**h3k27me3_bbdict},
                                        useSigma=False, all_ind_to_region=all_ind_to_region,
                                        ylimdict = {'H3K27me3 Treg' : (0, 5), 'H3K27ac Treg' : (0, 40)},
                                        bwcolor='black', gene_name_fontsize=14,
                                                                                    ax = ax)
        elif i1 < i2:
            ax = axs[c]
            _, ((ax1, ax2), bw_ax, place1, place2), ax = make_scoping_plot_for_fig1(treg_5kb, i1, i2, get_treg_mat, 
                                                                                    label, resolution=res, 
                                                                                    d_LR=d_LR,
                                                                                    d_LL=d_LL,
                                                                                    d_RR=d_RR,
                                                                                    d_RL=d_RL,
                                        ignore_set = all_gene_names_except_these, 
                                        vmin=vmin, vmax=vmax, 
                                        vert_bbdict={**h3k27me3_bbdict},
                                        useSigma=False, all_ind_to_region=all_ind_to_region,
                                        ylimdict = {'H3K27me3 Treg' : (0, 5), 'H3K27ac Treg' : (0, 40)},
                                        bwcolor='black', gene_name_fontsize=14,
                                                                                    ax = ax)

        elif i1 > i2:
            ax = axs[c]
            _, ((ax1, ax2), bw_ax, place1, place2), ax = make_scoping_plot_for_fig1(treg_5kb, i1, i2, get_tcon_mat, 
                                                                                    label, resolution=res, 
                                                                                    d_LR=d_LR,
                                                                                    d_LL=d_LL,
                                                                                    d_RR=d_RR,
                                                                                    d_RL=d_RL,
                                        ignore_set = all_gene_names_except_these, 
                                        vmin=vmin, vmax=vmax, 
                                        vert_bbdict={**h3k27me3_bbdict},
                                        useSigma=False, all_ind_to_region=all_ind_to_region,
                                        ylimdict = {'H3K27me3 Treg' : (0, 5), 'H3K27ac Treg' : (0, 40)},
                                        bwcolor='black', gene_name_fontsize=14,
                                                                                    ax = ax)
        ax.set_title("")
            

        if c != n*n-4:
            for a in bw_ax:
                a.set_yticklabels([])
            ax1.set_xlabel("")
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])    
        ax2.set_yticklabels([])
        ax2.spines['left'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)

        bw_ax[0].set_xticklabels([])
        ax2.set_xticklabels([])

        ax2.set_yticks([])
        ax1.set_xticks([])
        for a in bw_ax:
            # a.remove()
            pass
        
        if c % n != 0:
            ax1.remove()
        if n*n - c > n:
            for a in bw_ax:
                a.remove()
            ax2.remove()
            pass
        else:
            if c < n-1:
                for a in bw_ax:
                    a.set_ylabel("")
        if i1 == i2:
            lower_color = 'blue'
            upper_color = 'red'
        elif i2 > i1:
            lower_color = 'red'
            upper_color = 'red'
        elif i2 < i1:
            lower_color = 'blue'
            upper_color = 'blue'
        for _, spine in enumerate(axs[c].spines.keys()):
            if _ % 2 == 0:
                axs[c].spines[spine].set_edgecolor(lower_color)
            else:
                axs[c].spines[spine].set_edgecolor(upper_color)
            axs[c].spines[spine].set_linewidth(2)
            axs[c].spines[spine].set_visible(True)
        if i1 == i2:
            z = np.linspace(*axs[c].get_ylim(), 2000)
            axs[c].plot(z, z, color='white', linewidth=2)
            axs[c].plot(z[:-15], z[15:], color='red', linewidth=2)
            axs[c].plot(z[15:], z[:-15], color='blue', linewidth=2)
    return fig





def plot_active1_hub(all_ind_to_region, all_gene_names, get_treg_mat, get_tcon_mat, treg_5kb, tcon_5kb, deltas):
    all_gene_names_except_these = all_gene_names - set(['Pdlim7', 'Cd37', 'Cd27', 'Ikzf4', 'Lck', 'Jund', 'Cd4', 'Coro1b', 
            'Lta', 'Notch1', 'Id3', 'Cxcr5', 'Il27ra', 'Jak3', 'Lta', 'Nfkb2', 'Jak2', 'Tnfrsf10b', 'Notch1',
                                                        'Pdlim7', 'Cd37', 'Cd27', 'Ikzf4', 'Lck', 'Jund', 'Cd4', 'Coro1b', 
            'Lta', 'Notch1', 'Id3', 'Cxcr5', 'Il27ra', 'Jak3', 'Lta', 'Nfkb2', 'Jak2', 'Tnfrsf10b',
            'Irf3', 'Tnf', 'NCr3', 'Flt3l', 'Notch1']                                                    
                                                    )

    inds = [2670, 2696, 3887, 4166, 5264, 6099, 7299, 9010, 9629]
    deltas = [0, 100_000, 100_000, 0, 50_000, 50_000, -25_000, 0, 0, 0, 0]
    deltas = dict(zip(inds, deltas))

    d = 1_000_000
    res = 100_000
    n = len(inds)
    label = ''
    fig, axs = init_subplots_exact(n*n, n, fgsz = (2*30*mm, 2*30*mm), space=1.1)
    for c, (i1, i2) in enumerate(itertools.product(inds, inds)):
        if (i1 != i2) and (all_ind_to_region[i1][0] == all_ind_to_region[i2][0]):
            vmin = 1e-4
            vmax = 5e-3
        elif (i1 != i2):
            vmin = 1e-5
            vmax = 6e-4
        else:
            vmin = 1e-4
            vmax = 5e-1
        dL = deltas[i1]
        dR = deltas[i2]
        d_LR, d_LL, d_RR, d_RL = d-dL, d+dL, d-dR, d+dR
        if i1 == i2:
            ax = axs[c]
            _, ((ax1, ax2), bw_ax, place1, place2), ax = make_scoping_plot_for_fig1_near_diag(treg_5kb, tcon_5kb, i1, i2, get_treg_mat,
                                                                                        get_tcon_mat,
                                                                                    label, resolution=10_000, 
                                                                                    d_LR=d_LR,
                                                                                    d_LL=d_LL,
                                                                                    d_RR=d_RR,
                                                                                    d_RL=d_RL,
                                        ignore_set = all_gene_names_except_these, 
                                        vmin=vmin, vmax=vmax, 
                                        vert_bbdict={},
                                        useSigma=False, all_ind_to_region=all_ind_to_region,
                                        ylimdict = {'H3K27me3 Treg' : (0, 5), 'H3K27ac Treg' : (0, 40)},
                                        bwcolor='black', gene_name_fontsize=14,
                                                                                    ax = ax)
        elif i1 < i2:
            ax = axs[c]
            _, ((ax1, ax2), bw_ax, place1, place2), ax = make_scoping_plot_for_fig1(treg_5kb, i1, i2, get_treg_mat, 
                                                                                    label, resolution=res, 
                                                                                    d_LR=d_LR,
                                                                                    d_LL=d_LL,
                                                                                    d_RR=d_RR,
                                                                                    d_RL=d_RL,
                                        ignore_set = all_gene_names_except_these, 
                                        vmin=vmin, vmax=vmax, 
                                        vert_bbdict={},
                                        useSigma=False, all_ind_to_region=all_ind_to_region,
                                        ylimdict = {'H3K27me3 Treg' : (0, 5), 'H3K27ac Treg' : (0, 40)},
                                        bwcolor='black', gene_name_fontsize=14,
                                                                                    ax = ax)

        elif i1 > i2:
            ax = axs[c]
            _, ((ax1, ax2), bw_ax, place1, place2), ax = make_scoping_plot_for_fig1(treg_5kb, i1, i2, get_tcon_mat, 
                                                                                    label, resolution=res, 
                                                                                    d_LR=d_LR,
                                                                                    d_LL=d_LL,
                                                                                    d_RR=d_RR,
                                                                                    d_RL=d_RL,
                                        ignore_set = all_gene_names_except_these, 
                                        vmin=vmin, vmax=vmax, 
                                        vert_bbdict={},
                                        useSigma=False, all_ind_to_region=all_ind_to_region,
                                        ylimdict = {'H3K27me3 Treg' : (0, 5), 'H3K27ac Treg' : (0, 40)},
                                        bwcolor='black', gene_name_fontsize=14,
                                                                                    ax = ax)
            
        ax.set_title("")
            

        if c != n*n-4:
            for a in bw_ax:
                a.set_yticklabels([])
            ax1.set_xlabel("")
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])    
        ax2.set_yticklabels([])
        ax2.spines['left'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)

        bw_ax[0].set_xticklabels([])
        ax2.set_xticklabels([])

        ax2.set_yticks([])
        ax1.set_xticks([])
        for a in bw_ax:
            a.remove()
        
        if c % n != 0:
            ax1.remove()
        if n*n - c > n:
            # for a in bw_ax:
                # a.remove()
            ax2.remove()
        if i1 == i2:
            lower_color = 'blue'
            upper_color = 'red'
        elif i2 > i1:
            lower_color = 'red'
            upper_color = 'red'
        elif i2 < i1:
            lower_color = 'blue'
            upper_color = 'blue'
        for _, spine in enumerate(axs[c].spines.keys()):
            if _ % 2 == 0:
                axs[c].spines[spine].set_edgecolor(lower_color)
            else:
                axs[c].spines[spine].set_edgecolor(upper_color)
            axs[c].spines[spine].set_linewidth(2)
            axs[c].spines[spine].set_visible(True)
        if i1 == i2:
            z = np.linspace(*axs[c].get_ylim(), 2000)
            axs[c].plot(z, z, color='white', linewidth=2)
            axs[c].plot(z[:-15], z[15:], color='red', linewidth=2)
            axs[c].plot(z[15:], z[:-15], color='blue', linewidth=2)
    return fig