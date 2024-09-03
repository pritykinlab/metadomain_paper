from aux_functions import index_to_bedtool, get_pileup_from_bigwig
import numpy as np
from plotting_functions import init_subplots_exact
import matplotlib.pyplot as plt

def make_pileups(bed1, bed2, bwdict, delta=20_000, bins=2000,
                 label1='Invariant', label2='Variant', vmax = 100):
    protac_invariant_histone_dict = {}
    for key, bw in bwdict.items():
        pileup = get_pileup_from_bigwig(bw, index_to_bedtool(bed1), delta = delta, bins=bins)
        protac_invariant_histone_dict[key] = pileup

    protac_variant_histone_dict = {}
    for key, bw in bwdict.items():
        pileup = get_pileup_from_bigwig(bw, index_to_bedtool(bed2),  delta = delta, bins=bins)
        protac_variant_histone_dict[key] = pileup

    n = len(protac_invariant_histone_dict)
    fig, axs = init_subplots_exact(n, 1, dpi = 100)
    for c, (key, v) in enumerate(protac_invariant_histone_dict.items()):
        w = protac_variant_histone_dict[key]
        xs = np.linspace(-5e3, 5e3, w.shape[1])
        plt.sca(axs[c])
        plt.plot(xs, np.nanmean(v, axis=0), label=f'{label1}')
        plt.plot(xs, np.nanmean(w, axis=0), linestyle='--', label=f'{label2}')
        plt.title(key)
        plt.legend(handlelength=2)
        plt.xlabel("Distance from motif")


    for c, (key, v) in enumerate(protac_invariant_histone_dict.items()):
        w = protac_variant_histone_dict[key]
        fig, axs = init_subplots_exact(2, 1, dpi = 100)
        o1 = np.argsort(np.nanmean(-v[:, :], axis=1))
        o2 = np.argsort(np.nanmean(-w[:, :], axis=1))
        
        axs[0].matshow(v[o1, :], cmap='gist_heat_r', zorder=3, vmin=0, vmax=vmax, aspect='auto')
        axs[0].set_title(f'{key} at {label1} sites')

        axs[1].matshow(w[o2, :], cmap='gist_heat_r', zorder=3, vmin=0, vmax=vmax, aspect='auto')
        axs[1].set_title(f'{key} at {label2} sites')
        for ax in axs:
            plt.sca(ax)
            plt.xticks([])