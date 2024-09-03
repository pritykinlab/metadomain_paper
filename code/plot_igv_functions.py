from plotting_functions import *



def plot_igv_track(place, new_bwdict, histone_bwdict, peak_set, new_ylimdict, histone_ylimdict, label_for_title=None, alpha=.5,
                   outpath = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz',
                   linewidth=0, do_merge=True, chrom_converter = {}):
    n = len(new_bwdict) + len(histone_bwdict) + 1
    place = place
    fig, axs = init_subplots_exact(n, n, fgsz=(12, .5), sharex=True, dpi=100)
    ax = axs[0]
    plot_one_GTF_on_axis(place, ax, outpath=outpath)
    ax.set_yticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(bottom=False)
    for c, name in enumerate(new_bwdict):
        ax = axs[c + 1]
        plot_bigwig_on_axis_fill(add_chr_to_anc(place), 
                            new_bwdict, ax, name, ylimdict=new_ylimdict, label_for_title=name, alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))


    for c, name in enumerate(histone_bwdict):
        ax = axs[c + 1 + len(new_bwdict)]
        plot_bigwig_on_axis_fill(place, 
                            histone_bwdict, ax, name, ylimdict=histone_ylimdict, label_for_title=name,
                            color='blue', chrom_converter=chrom_converter,
                            alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))

    if peak_set is not None:    
        peaks = peak_set.intersect([add_chr_to_anc(place)], u=True).sort()
        if do_merge:
            peaks = peaks.merge()
        for ax in axs[1:]:
            seen = set()
            for peak in peaks:
                s, e = make_int(peak)[1:3]
                if (s, e) in seen:
                    continue
                else: 
                    seen.add((s, e))
                ax.axvspan(s, e, color='red', alpha=.2, linewidth=linewidth)
    for ax in axs:
        ax.grid(False)   
        ax.set_yticks([])
    for ax in axs[1:]:
        ax.set_xticks([])



def plot_igv_track_with_RNA(place, new_bwdict, histone_bwdict, rna_bwdict, peak_set, new_ylimdict, histone_ylimdict, label_for_title=None, alpha=.5,
                   outpath = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz',
                   linewidth=0, do_merge=True, chrom_converter = {}):
    n = len(new_bwdict) + len(histone_bwdict) + len(rna_bwdict) + 1
    place = place
    fig, axs = init_subplots_exact(n, n, fgsz=(12, .5), sharex=True, dpi=100, yspace=1.4)
    ax = axs[0]
    plot_one_GTF_on_axis(place, ax, outpath=outpath)
    ax.set_yticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(bottom=False)
    for c, name in enumerate(new_bwdict):
        ax = axs[c + 1]
        plot_bigwig_on_axis_fill(add_chr_to_anc(place), 
                            new_bwdict, ax, name, ylimdict=new_ylimdict, label_for_title=name, alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))


    for c, name in enumerate(histone_bwdict):
        ax = axs[c + 1 + len(new_bwdict)]
        plot_bigwig_on_axis_fill(place, 
                            histone_bwdict, ax, name, ylimdict=histone_ylimdict, label_for_title=name,
                            color='blue', chrom_converter=chrom_converter,
                            alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))

    if peak_set is not None:    
        peaks = peak_set.intersect([add_chr_to_anc(place)], u=True).sort()
        if do_merge:
            peaks = peaks.merge()
        for ax in axs[1:]:
            seen = set()
            for peak in peaks:
                s, e = make_int(peak)[1:3]
                if (s, e) in seen:
                    continue
                else: 
                    seen.add((s, e))
                ax.axvspan(s, e, color='red', alpha=.2, linewidth=linewidth)
    for ax in axs:
        ax.grid(False)   
        ax.set_yticks([])
    for ax in axs[1:]:
        ax.set_xticks([]) 
    print("HI")
    return axs



def plot_igv_track_with_RNA(place, new_bwdict, histone_bwdict, rna_bwdict, peak_set, new_ylimdict, histone_ylimdict, label_for_title=None, alpha=.5,
                   outpath = '/Genomics/argo/users/gdolsten/pritlab/jupys/tregs/pygbrowse/Mus_musculus.GRCm38.93.chr.gff3.gz',
                   linewidth=0, do_merge=True, chrom_converter = {}):
    n = len(new_bwdict) + len(histone_bwdict) + len(rna_bwdict) + 1
    place = place
    fig, axs = init_subplots_exact(n, n, fgsz=(12, .5), sharex=True, dpi=100, yspace=1.4)
    ax = axs[0]
    plot_one_GTF_on_axis(place, ax, outpath=outpath)
    ax.set_yticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(bottom=False)
    for c, name in enumerate(new_bwdict):
        ax = axs[c + 1]
        plot_bigwig_on_axis_fill(add_chr_to_anc(place), 
                            new_bwdict, ax, name, ylimdict=new_ylimdict, label_for_title=name, alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))


    for c, name in enumerate(histone_bwdict):
        ax = axs[c + 1 + len(new_bwdict)]
        plot_bigwig_on_axis_fill(place, 
                            histone_bwdict, ax, name, ylimdict=histone_ylimdict, label_for_title=name,
                            color='blue', chrom_converter=chrom_converter,
                            alpha=1)

        for ax in axs:
            ax.ticklabel_format(axis='x', style='sci', scilimits=(6, 6))

    if peak_set is not None:    
        peaks = peak_set.intersect([add_chr_to_anc(place)], u=True).sort()
        if do_merge:
            peaks = peaks.merge()
        for ax in axs[1:]:
            seen = set()
            for peak in peaks:
                s, e = make_int(peak)[1:3]
                if (s, e) in seen:
                    continue
                else: 
                    seen.add((s, e))
                ax.axvspan(s, e, color='red', alpha=.2, linewidth=linewidth)
    for ax in axs:
        ax.grid(False)   
        ax.set_yticks([])
    for ax in axs[1:]:
        ax.set_xticks([]) 
    print("HI")
    return axs