from coolbox.api import ChromName, XAxis, Title, Spacer, BED, HiCMat, BEDPE, Pairs, ArcsCoverage
import numpy as np

def save_plot_with_coolbox_diag(ind, d1, d2, all_ind_to_region, cools_frame, end='svg', **kwargs):
    fig, L, R = plot_diag_with_coolbox(ind, d1, d2, all_ind_to_region, frame2=cools_frame, **kwargs)
    return fig, L, R

def plot_diag_with_coolbox(ind, d1, d2, all_ind_to_region, treg_coolpath='./Treg_all.mcool', 
                               tcon_coolpath='./Treg_all.mcool', 
                               outdir='./plots/', frame2=None, 
                               vmin = 5e-4, vmax=5e-1, depth_ratio=.45,
                               ctcf_paths = ['./peaks/motifs/all_ctcf_peaks.txt'],
                          print_p = False, resolution=250_000, 
                          xaxis_width = 20, cool_width = 25, hic_height = 8, 
                          treg_label = 'Treg Hi-C', tcon_label = 'Tcon Hi-C'):

    cool_1 = (ChromName() + XAxis(width=xaxis_width) 
                 + HiCMat(treg_coolpath, balance=True, cmap='gist_heat_r', 
                    resolution=5_000, depth_ratio=depth_ratio, 
                    min_value = np.log10(vmin), 
                          max_value = np.log10(vmax),
                           transform = 'log10', height=hic_height,
                )  + Title(treg_label)

             + (Pairs(f'./arcs_for_coolbox/thresh=0_ns_loops.pairs', line_width=1.5, orientation='inverted', color='lightgray', height=1)
                + ArcsCoverage(f'./arcs_for_coolbox/thresh=0_treg_loops.pairs', line_width=1.5, orientation='inverted', color='red')
               ) + Title("Loops")
             + BED('./FINAL_boundaries/all_boundaries.bed', display='collapsed', border_width=0, border_color='black', color='black', properties_dict={"alpha" : .2}, height=.13, labels=False, show_score=False) + Title("Boundaries")
             + Spacer(.25) 

            )
    cool_2 = (
             HiCMat(tcon_coolpath, balance=True, cmap='gist_heat_r', 
                resolution=5_000, depth_ratio=depth_ratio, min_value = np.log10(vmin), 
                      max_value = np.log10(vmax), transform = 'log10', height=hic_height,
                )   + Title(tcon_label)
             + (Pairs(f'./arcs_for_coolbox/thresh=0_ns_loops.pairs', line_width=1.5, orientation='inverted', color='lightgray', height=1)
                + ArcsCoverage(f'./arcs_for_coolbox/thresh=0_tcon_loops.pairs', line_width=1.5, orientation='inverted', color='blue')
               ) + Title("Loops")
             + BED('./FINAL_boundaries/all_boundaries.bed', display='collapsed', border_width=0, border_color='black', color='black', properties_dict={"alpha" : .2}, height=.13, labels=False, show_score=False) + Title("Boundaries")
             + Spacer(.25) 
            )
    
    for name, ctcf_path in ctcf_paths.items():
        cool_1 += (
                 BED(ctcf_path, border_color='black', 
                          color='none', display='collapsed', height=0.3, 
                           labels=False, show_score=False) + Title(f"CTCF Motif {name}")
				+ Spacer(.25)
                ) 

        cool_2 += (
                BED(ctcf_path, border_color='black', 
                          color='none', display='collapsed', height=0.3, 
                           labels=False, show_score=False) + Title(f"CTCF Motif {name}")
				+ Spacer(.25)
                )

    cool = cool_1 + cool_2
    cool += frame2
    chrom1, s1, e1 = all_ind_to_region[ind]
    s1 = int(s1 - d1*resolution)
    e1 = int(e1 + d2*resolution)
    cool.properties['width'] = cool_width
    fig = cool.plot(f"chr{chrom1}:{s1}-{e1}")
    return fig, (chrom1, s1, e1), (chrom1, s1, e1)




def plot_diag_with_coolbox_region(reg1, treg_coolpath='./Treg_all.mcool', 
                               tcon_coolpath='./Treg_all.mcool', 
                               outdir='./plots/', frame2=None, 
                               vmin = 5e-4, vmax=5e-1, depth_ratio=.45,
                               ctcf_paths = ['./peaks/motifs/all_ctcf_peaks.txt'],
                          print_p = False, resolution=250_000, 
                          xaxis_width = 20, cool_width = 25, hic_height = 8):

    cool_1 = (ChromName() + XAxis(width=xaxis_width) 
                 + HiCMat(treg_coolpath, balance=True, cmap='gist_heat_r', 
                    resolution=5_000, depth_ratio=depth_ratio, min_value = np.log10(vmin), 
                          max_value = np.log10(vmax), transform = 'log10', height=hic_height,
                )  + Title("Treg Hi-C")

             + (Pairs(f'./arcs_for_coolbox/thresh=0_ns_loops.pairs', line_width=1.5, orientation='inverted', color='lightgray', height=1)
                + ArcsCoverage(f'./arcs_for_coolbox/thresh=0_treg_loops.pairs', line_width=1.5, orientation='inverted', color='red')
               ) + Title("Loops")
             + BED('./FINAL_boundaries/all_boundaries.bed', display='collapsed', color='black', properties_dict={"alpha" : .2}, height=.13, labels=False, show_score=False) + Title("Boundaries")
             + Spacer(.25) 

            )
    cool_2 = (
             HiCMat(tcon_coolpath, balance=True, cmap='gist_heat_r', 
                resolution=5_000, depth_ratio=depth_ratio, min_value = np.log10(vmin), 
                      max_value = np.log10(vmax), transform = 'log10', height=hic_height,
                )   + Title("Tcon Hi-C")
             + (Pairs(f'./arcs_for_coolbox/thresh=0_ns_loops.pairs', line_width=1.5, orientation='inverted', color='lightgray', height=1)
                + ArcsCoverage(f'./arcs_for_coolbox/thresh=0_tcon_loops.pairs', line_width=1.5, orientation='inverted', color='blue')
               ) + Title("Loops")
             + BED('./FINAL_boundaries/all_boundaries.bed', display='collapsed', color='black', properties_dict={"alpha" : .2}, height=.13, labels=False, show_score=False) + Title("Boundaries")
             + Spacer(.25) 
            )
    
    for name, ctcf_path in ctcf_paths.items():
        cool_1 += (
                 BED(ctcf_path, border_color='black', 
                          color='none', display='collapsed', height=0.3, 
                           labels=False, show_score=False) + Title(f"CTCF Motif {name}")
				+ Spacer(.25)
                ) 

        cool_2 += (
                BED(ctcf_path, border_color='black', 
                          color='none', display='collapsed', height=0.3, 
                           labels=False, show_score=False) + Title(f"CTCF Motif {name}")
				+ Spacer(.25)
                )

    cool = cool_1 + cool_2
    cool += frame2
    cool.properties['width'] = cool_width
    fig = cool.plot(f"{reg1[0]}:{reg1[1]}-{reg1[2]}")
    return fig, reg1, reg1


