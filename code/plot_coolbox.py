from coolbox.api import ChromName, XAxis, Title, Spacer, BED, HiCMat
import numpy as np

def save_plot_with_coolbox_diag(i1, i2, d1, d2, cools_frame, **kwargs):
    ind1, ind2 = (i1-d1, i1+d1), (i2-d2, i2+d2+1)
    assert i1 == i2
    fig, L, R = plot_diag_with_coolbox(ind1, ind2, frame2=cools_frame, **kwargs)
#     fig.savefig(f"./plots/fig2_all_intra/diag_{i1}_h3k27me3_tcon.{end}")
    return fig, L, R

def plot_diag_with_coolbox(ind1, ind2, all_ind_to_region, treg_coolpath='./Treg_all.mcool', 
                               tcon_coolpath='./Treg_all.mcool', 
                               outdir='./plots/', frame2=None, 
                               vmin = 5e-4, vmax=5e-1, depth_ratio=.45,
                               ctcf_paths = ['./peaks/motifs/all_ctcf_peaks.txt'],
                          print_p = False):

    cool_1 = (ChromName() + XAxis() 
                 + HiCMat(treg_coolpath, balance=True, cmap='gist_heat_r', 
                    resolution=5_000, depth_ratio=depth_ratio, min_value = np.log10(vmin), 
                          max_value = np.log10(vmax), transform = 'log10',
                )  + Title("Treg Hi-C")
             + Spacer(.25) 
            )
    cool_2 = (
             HiCMat(tcon_coolpath, balance=True, cmap='gist_heat_r', 
                resolution=5_000, depth_ratio=depth_ratio, min_value = np.log10(vmin), 
                      max_value = np.log10(vmax), transform = 'log10',
                )   + Title("Tcon Hi-C")
             + Spacer(.25) 
            )
    
    for name, ctcf_path in ctcf_paths.items():
        cool_1 += (
                 BED(ctcf_path, border_color='black', 
                          color=None, display='collapsed', height=0.1, 
                           labels=False, show_score=False) + Title(f"CTCF Motif {name}")
                )

        cool_2 += (
                BED(ctcf_path, border_color='black', 
                          color=None, display='collapsed', height=0.1, 
                           labels=False, show_score=False) + Title(f"CTCF Motif {name}")
                )

    cool = cool_1 + cool_2
    cool += frame2
    chrom1, s1, _ = all_ind_to_region[ind1[0]]
    chrom1, _, e1 = all_ind_to_region[ind1[1]]

    chrom2, s2, _ = all_ind_to_region[ind2[0]]
    chrom2, _, e2 = all_ind_to_region[ind2[1]]
    
    fig = cool.plot(f"chr{chrom1}:{s1}-{e2}")
    return fig, (chrom1, s1, e1), (chrom2, s2, e2)