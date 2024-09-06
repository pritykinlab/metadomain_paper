import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd

import seaborn as sns

def make_decay_plot(expected_dict):
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    axs = np.ravel(axs)
    for i, condition in enumerate(list(expected_dict.keys()) + ['All']):
        ax = axs[i]
        if i == 0:
            c1 = condition
        if i == 1:
            c2 = condition
        if i != 2:
            df = expected_dict[condition]
            for chrom, x in df.groupby('region'):
                if chrom == 'M' or chrom=='Y':
                    continue
                ax.plot(x.diag, x['balanced.avg'], label=chrom)
            ax.set_yscale('log')
            ax.set_xscale('log')
            # ax.legend()
            ax.set_title(condition)
            ax.set_ylim([1e-10, 5e-2])

        if i == 2:
            df1 = expected_dict[c1]
            df2 = expected_dict[c2]
            for chrom in np.unique(df.region.values):
                x = df1[df1.region == chrom].sort_values('diag')['balanced.avg'].values
                y = df2[df2.region == chrom].sort_values('diag')['balanced.avg'].values
                ax.plot(np.log2(x)-np.log2(y))
            # ax.set_yscale('log')
            ax.set_xscale('log')

    fig.savefig('./plots/qc/decay/decay_balanced.png')

        
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    axs = np.ravel(axs)
    for i, condition in enumerate(list(expected_dict.keys()) + ['All']):
        ax = axs[i]
        if i == 0:
            c1 = condition
        if i == 1:
            c2 = condition
        if i != 2:
            df = expected_dict[condition]
            for chrom, x in df.groupby('region'):
                if chrom == 'M' or chrom=='Y':
                    continue
                ax.plot(x.diag, x['count.avg'], label=chrom)
            ax.set_yscale('log')
            ax.set_xscale('log')
            # ax.legend()
            ax.set_title(condition)
            ax.set_ylim([1e-4, 1e3])

        if i == 2:
            df1 = expected_dict[c1]
            df2 = expected_dict[c2]
            for chrom in np.unique(df.region.values):
                x = df1[df1.region == chrom].sort_values('diag')['count.avg'].values
                y = df2[df2.region == chrom].sort_values('diag')['count.avg'].values
                ax.plot(np.log2(x)-np.log2(y))
            # ax.set_yscale('log')
            ax.set_xscale('log')
        # ax.legend()
        ax.set_title(condition)
    fig.savefig('./plots/qc/decay/decay_count.png')



def atac_seq_bars(dist_df, file):
    df = pd.DataFrame()
    names = [] 
    vs = []
    xs = []
    for name in ['Non-diff loops']:
        tr, tl, br, bl, = dist_df[name]
        names += [name]*4

        total = np.sum([tr, tl, br, bl])
        v = [tr/total, tl/total, br/total, bl/total]
        vs = vs + v
        xs = xs + ['++', '+-', '-+', '- -']
    for name in ['Treg loops']:
        tr, tl, br, bl, = dist_df[name]
        names += [name]*4

        total = np.sum([tr, tl, br, bl])
        v = [tr/total, tl/total, br/total, bl/total]
        vs = vs + v
        xs = xs + ['++', '+-', '-+', '- -']
    
    for name in ['Tcon loops']:
        tr, tl, br, bl, = dist_df[name]
        names += [name]*4

        total = np.sum([tr, tl, br, bl])
        v = [tr/total, tl/total, br/total, bl/total]
        vs = vs + v
        xs = xs + ['++', '+-', '-+', '- -']

    # print(xs, vs, names)
    df['xs'] = xs
    df['vs'] = vs
    df['name'] = names
    print(df.name)
    fig, ax  = plt.subplots()
    sns.barplot(x='xs', y='vs', hue='name', data=df, ax=ax)
    ax.set_xlabel("ATAC FC in both anchors: Tcon ÷ Treg")
    ax.set_ylabel("Percentage of loops")
    fig.savefig(f'./plots/atac_bar/{file}')


def atac_seq_cdfs(dist_df, file):
    df = pd.DataFrame()
    names = [] 
    vs = []
    xs = []
    for name in ['mustache_tcon_loops']:
        tr, tl, br, bl, = dist_df[name]
        names += [name]*4

        total = np.sum([tr, tl, br, bl])
        v = [tr/total, tl/total, br/total, bl/total]
        vs = vs + v
        xs = xs + ['++', '+-', '-+', '- -']
    for name in ['mustache_nonsig_loops']:
        tr, tl, br, bl, = dist_df[name]
        names += [name]*4

        total = np.sum([tr, tl, br, bl])
        v = [tr/total, tl/total, br/total, bl/total]
        vs = vs + v
        xs = xs + ['++', '+-', '-+', '- -']
    
    for name in ['mustache_treg_loops']:
        tr, tl, br, bl, = dist_df[name]
        names += [name]*4

        total = np.sum([tr, tl, br, bl])
        v = [tr/total, tl/total, br/total, bl/total]
        vs = vs + v
        xs = xs + ['++', '+-', '-+', '- -']

    # print(xs, vs, names)
    df['xs'] = xs
    df['vs'] = vs
    df['name'] = names
    print(df.name)
    fig, ax  = plt.subplots()
    sns.barplot(x='xs', y='vs', hue='name', data=df, ax=ax)
    ax.set_xlabel("ATAC FC in both anchors: Tcon ÷ Treg")
    ax.set_ylabel("Percentage of loops")
    fig.savefig(f'./plots/atac_bar/{file}')    



def atac_figures(xs_tcon, ys_tcon, xs_kiko, ys_kiko, cs, name):
    fig, axs = plt.subplots(1, 3, figsize=(21, 5))
    xs = xs_tcon
    ys = ys_tcon

    ax = axs[0]
    c = ax.scatter(xs, ys, c=cs, s= 4, vmin=-2, vmax=2, cmap=cm.bwr)
    # plt.colorbar(c)
    ax.set_xlabel("ATAC signal, Anchor 1: Tcon ÷ Treg")
    ax.set_ylabel("ATAC signal, Anchor 2: Tcon ÷ Treg")
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3, 3])
    ax.vlines(0, -3, 3, linestyle='--', color='gray', alpha=.5)
    ax.hlines(0, -3, 3, linestyle='--', color='gray', alpha=.5)

    topleft = (xs < 0) & (ys > 0)
    topright = (xs > 0) & (ys > 0)
    bottomleft = (xs < 0) & (ys < 0)
    bottomright = (xs > 0) & (ys < 0)

    from copy import deepcopy
    poop = deepcopy([topright.sum(), topleft.sum(), bottomright.sum(), bottomleft.sum(), ])

    ax.text(0.1, 0.9, f'{topleft.sum()}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.9, 0.9, f'{topright.sum()}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.1, 0.1, f'{bottomleft.sum()}', size=10, color='purple', transform=ax.transAxes)    
    ax.text(0.9, 0.1, f'{bottomright.sum()}', size=10, color='purple', transform=ax.transAxes)

    ax.text(0.1, 0.85, f'{np.round(topleft.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.9, 0.85, f'{np.round(topright.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.1, 0.05, f'{np.round(bottomleft.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)    
    ax.text(0.9, 0.05, f'{np.round(bottomright.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)


    xs = xs_kiko
    ys = ys_kiko

    ax = axs[1]
    c = ax.scatter(xs, ys, c=cs, s= 4, vmin=-2, vmax=2, cmap=cm.bwr)
    # plt.colorbar(c)
    ax.set_xlabel("ATAC signal, Anchor 1: KIKO ÷ Treg")
    ax.set_ylabel("ATAC signal, Anchor 2: KIKO ÷ Treg")
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3, 3])
    ax.vlines(0, -3, 3, linestyle='--', color='gray', alpha=.5)
    ax.hlines(0, -3, 3, linestyle='--', color='gray', alpha=.5)

    topleft = (xs < 0) & (ys > 0)
    topright = (xs > 0) & (ys > 0)
    bottomleft = (xs < 0) & (ys < 0)
    bottomright = (xs > 0) & (ys < 0)
    # poop = deepcopy([topright.sum(), topleft.sum(), bottomright.sum(), bottomleft.sum(), ])


    ax.text(0.1, 0.9, f'{topleft.sum()}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.9, 0.9, f'{topright.sum()}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.1, 0.1, f'{bottomleft.sum()}', size=10, color='purple', transform=ax.transAxes)    
    ax.text(0.9, 0.1, f'{bottomright.sum()}', size=10, color='purple', transform=ax.transAxes)

    ax.text(0.1, 0.85, f'{np.round(topleft.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.9, 0.85, f'{np.round(topright.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.1, 0.05, f'{np.round(bottomleft.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)    
    ax.text(0.9, 0.05, f'{np.round(bottomright.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)



    xs = xs_kiko - xs_tcon
    ys = ys_kiko - ys_tcon
    ax = axs[2]
    c = ax.scatter(xs, ys, c=cs, s= 4, vmin=-2, vmax=2, cmap=cm.bwr)
    # plt.colorbar(c)
    ax.set_xlabel("ATAC signal, Anchor 1: KIKO ÷ Tcon")
    ax.set_ylabel("ATAC signal, Anchor 2: KIKO ÷ Tcon")
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3, 3])
    ax.vlines(0, -3, 3, linestyle='--', color='gray', alpha=.5)
    ax.hlines(0, -3, 3, linestyle='--', color='gray', alpha=.5)

    topleft = (xs < 0) & (ys > 0)
    topright = (xs > 0) & (ys > 0)
    bottomleft = (xs < 0) & (ys < 0)
    bottomright = (xs > 0) & (ys < 0)

    ax.text(0.1, 0.9, f'{topleft.sum()}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.9, 0.9, f'{topright.sum()}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.1, 0.1, f'{bottomleft.sum()}', size=10, color='purple', transform=ax.transAxes)    
    ax.text(0.9, 0.1, f'{bottomright.sum()}', size=10, color='purple', transform=ax.transAxes)

    ax.text(0.1, 0.85, f'{np.round(topleft.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.9, 0.85, f'{np.round(topright.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)
    ax.text(0.1, 0.05, f'{np.round(bottomleft.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)    
    ax.text(0.9, 0.05, f'{np.round(bottomright.sum()/len(topleft),2)}', size=10, color='purple', transform=ax.transAxes)
    # poop = deepcopy([topright.sum(), topleft.sum(), bottomright.sum(), bottomleft.sum(), ])    
    fig.savefig(f'./plots/atac_scatterplot/{name}')
    return poop


# import cooltools.lib.numutils as nu 
def loop_pileups(tfs, loopsets, treg_expecteds, tcon_expecteds, cool_treg, cool_tconv):
    res=5000
    padfrac = 2
    import os
    for loopname in loopsets:

        for tf_name in tfs:
            os.makedirs(f'./plots/pileup_loops/', exist_ok=True)
            loopset = loopsets[loopname]
            tf = tfs[tf_name]
            seen = {}
            fcs = []
            vals_reg = []
            vals_cool = []
            for i in loopset.pair_to_bed(tf.slop(b=0, g='./annotations/chromsizes'), type='either'):
                place = tuple(i[:6])
                if seen.get(place, 0):
                    continue
                seen[place] = 1
                chrom, s = i[:2]
                e = i[5]

                s, e = map(int, [s, e])
                delta = e-s

                if (e - s)/res < 20:
                    continue
                try:
                    v_reg = cool_treg.matrix(balance=False).fetch((chrom, s-delta*padfrac, e+delta*padfrac))
                    v_cool = cool_tconv.matrix(balance=False).fetch((chrom, s-delta*padfrac, e+delta*padfrac))
                except Exception as e:
                    if 'Genomic region out of bounds' in str(e):
                        continue
                    else:
                        raise Exception
                
                n = v_reg.shape[0]
                
             
                treg_fc = np.log2((v_reg+1)/(treg_expecteds[chrom][:n, :n]+1))
                tcon_fc = np.log2((v_cool+1)/(tcon_expecteds[chrom][:n, :n]+1))
                
                treg_fc = nu.zoom_array(treg_fc, (100, 100))
                tcon_fc = nu.zoom_array(tcon_fc, (100, 100))
                fc = treg_fc - tcon_fc

                vals_reg.append(treg_fc)
                vals_cool.append(tcon_fc)

            fig, ax = plt.subplots(figsize=(8, 8))
            c = ax.matshow(np.nanmean(vals_reg, axis=0), cmap=cm.bwr, vmin=-1, vmax=1)        
            plt.colorbar(c)
            fig.savefig(f'plots/pileup_loops/{tf_name}_{loopname}_treg.png')


            fig, ax = plt.subplots(figsize=(8, 8))
            c = ax.matshow(np.nanmean(vals_cool, axis=0), cmap=cm.bwr, vmin=-1, vmax=1)        
            plt.colorbar(c)
            fig.savefig(f'plots/pileup_loops/{tf_name}_{loopname}_tcon.png')            



# import cooltools.lib.numutils as nu 
def superenhancer_pileups(tfs, treg_expecteds, tcon_expecteds, cool_treg, cool_tconv):
    res=5000
    padfrac = 4
    import os
    for tf_name in tfs:
        os.makedirs(f'./plots/loop_pileups/{tf_name}', exist_ok=True)
        tf = tfs[tf_name]
        seen = {}
        fcs = []
        vals_reg = []
        vals_cool = []
        
        for i in tf:
            place = tuple(i[:3])
            if seen.get(place, 0):
                continue
            seen[place] = 1
            chrom, s = i[:2]
            e = i[2]

            s, e = map(int, [s, e])
            delta = e-s

            if (e - s)/res < 20:
                continue
            try:
                v_reg = cool_treg.matrix(balance=False).fetch((chrom, s-delta*padfrac, e+delta*padfrac))
                v_cool = cool_tconv.matrix(balance=False).fetch((chrom, s-delta*padfrac, e+delta*padfrac))
            except Exception as e:
                if 'Genomic region out of bounds' in str(e):
                    continue
                else:
                    raise Exception
            n = v_reg.shape[0]
            treg_fc = (v_reg+1)/(treg_expecteds[chrom][:n, :n]+1)
            tcon_fc = (v_cool+1)/(tcon_expecteds[chrom][:n, :n]+1)            
            
            treg_fc = nu.zoom_array(treg_fc, (500, 500))
            tcon_fc = nu.zoom_array(tcon_fc, (500, 500))
            fc = treg_fc / tcon_fc

            vals_reg.append(treg_fc)
            vals_cool.append(tcon_fc)
            fcs.append(fc)
            
        fig, ax = plt.subplots(figsize=(8, 8))
        c = ax.matshow(np.nanmean(np.log2(fcs), axis=0), cmap=cm.bwr, vmin=-1, vmax=1)  
        ax.set_title(f"{fcs.shape[0]}")
        plt.colorbar(c)
        fig.savefig(f'plots/superenhancer_pileups/fc.png')


        fig, ax = plt.subplots(figsize=(8, 8))
        c = ax.matshow(np.nanmean(np.log2(vals_reg), axis=0), cmap=cm.bwr, vmin=-1, vmax=1)        
        ax.set_title(f"{vals_reg.shape[0]}")
        plt.colorbar(c)
        fig.savefig(f'plots/superenhancer_pileups/treg.png')


        fig, ax = plt.subplots(figsize=(8, 8))
        c = ax.matshow(np.nanmean(np.log2(vals_cool), axis=0), cmap=cm.bwr, vmin=-1, vmax=1)        
        ax.set_title(f"{vals_cool.shape[0]}")
        plt.colorbar(c)
        fig.savefig(f'plots/superenhancer_pileups/tcon.png')                    

def insulation_loop_heatmap(anchors_of_interest_dict, tfs, treg_insulation, tcon_insulation, chromlist):
    v_dict = {}
    for anc_name in anchors_of_interest_dict.keys():
        loops_of_interest = anchors_of_interest_dict[anc_name] 
        for tf_name in tfs:
            vs_plus = []
            vs_minus = []
            vs_both = []
            vs = []
            
            tf, sub = tfs[tf_name]
            seen = {}
            to_consider = list(loops_of_interest.pair_to_bed(tf.slop(b=0, g='./annotations/chromsizes'), type='both'))
            places = []
            for _, i in enumerate(to_consider):
                place = tuple(i[:6])
                if seen.get(place):
                    continue
                else:
                    seen[place] = 1
                l1, l2 = place[:3], place[3:6]
                chrom, s, e = l1[0], l1[1], l2[1]
                if chrom not in chromlist:
                    continue

                if s == e:
                    continue
                s, e = map(int, [s, e])
                delta = e-s
                L, R = s-delta*.5, e+delta*.5
                places.append((chrom, L, R))

            chroms, starts, ends = list(zip(*places))
            vs_treg = treg_insulation.stackup(chroms, starts, ends, missing=np.nan, oob=np.nan, bins=1000)
            v_dict[tf_name] = vs_treg

            chroms, starts, ends = list(zip(*places))
            vs_tcon = tcon_insulation.stackup(chroms, starts, ends, missing=np.nan, oob=np.nan, bins=1000)
            
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            axs[0].plot(np.nanmean(vs_treg, axis=0), label=f'treg: {vs_treg.shape[0]}')
            axs[1].plot(np.nanmean(vs_tcon, axis=0), label=f'tcon: {vs_tcon.shape[0]}')
            axs[0].set_ylim([-.1,1.3])
            axs[1].set_ylim([-.1,1.3])
            axs[0].legend()
            axs[1].legend()
            fig.suptitle(f'{anc_name}_{tf_name}.png')
            plt.tight_layout()
            fig.savefig(f'./plots/insulation_loops_heatmap/{tf_name}.png')
            

            vs_treg[np.isnan(vs_treg)] = 0        
            heatmap = sns.clustermap(vs_treg, col_cluster=False, vmin=-.6, vmax = .6, cmap=cm.bwr, )
            
            order_row = deepcopy(heatmap.dendrogram_row.reordered_ind)
            heatmap.savefig(f'./plots/insulation_loops_heatmap/{tf_name}_treg.png')        
            
            vs_tcon[np.isnan(vs_tcon)] = 0
            heatmap = sns.clustermap(vs_tcon[order_row, :], col_cluster=False, row_cluster=False, vmin=-.6, vmax = .6, cmap=cm.bwr, )
            heatmap.savefig(f'./plots/insulation_loops_heatmap/{tf_name}_tcon.png')        



def insulation_tf_heatmap(mustache_full_anchors, tfs, treg_insulation, tcon_insulation, chromlist):
    for tf_name in tfs:
        places = []

        tf = tfs[tf_name][0]
        sub = tfs[tf_name][1]
        for i in tf.slop(b=5000, g='./annotations/chromsizes'
                        ).subtract(sub.slop(b=5000, g='./annotations/chromsizes'), A=True).intersect(mustache_full_anchors, u=True):
            place = tuple(i[:3])
            chrom, s, e = place
            if chrom not in chromlist:
                continue
            s, e = map(int, [s, e])
            s -= 5000*40
            e += 5000*40
            places.append((chrom, s, e))
        chroms, starts, ends = list(zip(*places))
        vs_plus_treg = treg_insulation.stackup(chroms, starts, ends, missing=np.nan, oob=np.nan, bins=500)
        vs_plus_tcon = tcon_insulation.stackup(chroms, starts, ends, missing=np.nan, oob=np.nan, bins=500)

        fig, axs = plt.subplots(1, 2, figsize=(8, 4))
        axs[0].plot(np.nanmean(vs_plus_treg, axis=0))
        axs[1].plot(np.nanmean(vs_plus_tcon, axis=0))    
        fig.savefig(f'./plots/insulation_heatmap_tf/{tf_name}_lineplot.png')    

        v = vs_plus_treg
        v[np.isnan(v)] = 0
        heatmap = sns.clustermap(vs_plus_treg, col_cluster=False, vmin=-1.2, vmax = 1.2, cmap=cm.bwr, )
        order_row = deepcopy(heatmap.dendrogram_row.reordered_ind)
        heatmap.savefig(f'./plots/insulation_heatmap_tf/{tf_name}_treg_insulation.png')        

        v = vs_plus_tcon
        v[np.isnan(v)] = 0
        heatmap = sns.clustermap(vs_plus_tcon[order_row, :], col_cluster=False, row_cluster=False, vmin=-1.2, vmax = 1.2, cmap=cm.bwr, )
        heatmap.savefig(f'./plots/insulation_heatmap_tf/{tf_name}_tcon_insulation.png')        

import sklearn
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns
from copy import deepcopy
def pileups_final(tfs, cool_treg, mustache_full_anchors, chromlist, anchor_direction, treg_expecteds_balanced, treg_expecteds):
    n=100  
    rawdict_balanced = {}
    rawdict = {}
    chroms = {}
    for tf_name in tfs:
        print(tf_name)
        lists_of_chroms = []
        tf, sub = tfs[tf_name]
        raw = []
        raw_balanced = []
        for i in mustache_full_anchors.intersect(tf, u=True).subtract(sub, A=True)[:]:
            chrom, s, e = i[:3]
            if chrom not in chromlist:
                print(chrom, chromlist)
                continue
            s, e = map(int, [s, e])

            try:
                mat = np.diag((cool_treg.matrix(balance=True).fetch((chrom, s-5000*n, e+5000*n)))[::-1])
            except Exception as e:
                if 'Genomic region out of bounds' in str(e):
                    continue
                else:
                    print(e)
                    raise Exception  
            rawmat = (cool_treg.matrix(balance=False).fetch((chrom, s-5000*n, e+5000*n)))
            if rawmat.shape != (203, 203):
                print(chrom, s, e)
                continue
            balmat = (cool_treg.matrix(balance=True).fetch((chrom, s-5000*n, e+5000*n)))    
            assert balmat.shape == (203, 203)
            if '-' in anchor_direction[(chrom, str(s), str(e))]:
                rawmat = np.rot90(rawmat, k=2)
                balmat = np.rot90(balmat, k=2)
            raw.append(rawmat)
            raw_balanced.append(balmat)
            lists_of_chroms.append(chrom)
        chroms[tf_name] = lists_of_chroms
        rawdict[tf_name] = raw
        rawdict_balanced[tf_name] = raw_balanced

    vals_raw = {}
    vals_balanced = {}
    vals_diag_balanced_25 = {}
    vals_diag_balanced_50 = {}
    vals_diag_raw_25 = {}
    vals_diag_raw_50 = {}

    pc_rat = .5
    dict_types = {
        'balanced' : rawdict_balanced,
        'raw' : rawdict
    }

    mats_raw = {}
    mats_balanced = {}
    for t in dict_types:
        if t == 'balanced':
            for key in rawdict_balanced:  
                vals_balanced[key] = []
                vals_diag_balanced_25[key]  = []
                vals_diag_balanced_50[key]  = []

                mats_balanced[key] = []
                for _, i in enumerate(rawdict_balanced[key]):
                    expected = treg_expecteds_balanced[chroms[key][_]][:2*n+3, :2*n+3]
                    mat = np.log2((i+expected*pc_rat)/(expected+expected*pc_rat))
                    mats_balanced[key].append(mat)

                    crossdiag = np.diag(mat[::-1])
                    diag25 = np.diag(mat, k=25)
                    diag50 = np.diag(mat, k = 50)
                    vals_balanced[key].append(crossdiag)
                    vals_diag_balanced_25[key].append(diag25)
                    vals_diag_balanced_50[key].append(diag25)
                vals_balanced[key] = np.asarray(vals_balanced[key])
                vals_diag_balanced_25[key] = np.asarray(vals_diag_balanced_25[key])
                vals_diag_balanced_50[key] = np.asarray(vals_diag_balanced_50[key])
        if t == 'raw':
            for key in rawdict:  
                vals_raw[key] = []
                vals_diag_raw_25[key]  = []   
                vals_diag_raw_50[key]  = []               
                mats_raw[key] = []
                
                for _, i in enumerate(rawdict[key]):
                    expected = treg_expecteds[chroms[key][_]][:2*n+3, :2*n+3]
                    mat = np.log2((i+expected*pc_rat)/(expected+expected*pc_rat))
                    mats_raw[key].append(mat)
                    crossdiag = np.diag(mat[::-1])
                    diag25 = np.diag(mat, k = 25)
                    diag50 = np.diag(mat, k = 50)
                    vals_raw[key].append(crossdiag)
                    vals_diag_raw_25[key].append(diag25)
                    vals_diag_raw_50[key].append(diag50)

                vals_raw[key] = np.asarray(vals_raw[key])
                vals_diag_raw_25[key] = np.asarray(vals_diag_raw_25[key])
                vals_diag_raw_50[key] = np.asarray(vals_diag_raw_50[key])    
    c = 0
    for tf_name in vals_diag_raw_50:
        mat_raw = deepcopy(vals_diag_raw_50[tf_name])
        mat_balanced = deepcopy(vals_diag_balanced_50[tf_name])

        mat_raw[np.isnan(mat_raw)] = 0
        mat_balanced[np.isnan(mat_balanced)] = 0
        
        
        agg = AgglomerativeClustering(n_clusters = 6)
        agg = agg.fit((mat_raw))
        order = np.argsort(agg.labels_)

        vs = []
        for i in range(np.max(agg.labels_)+1):
            vs.append(np.mean(mat_raw[:, 80:120][agg.labels_ == i]))
        sortdict = dict(zip(np.argsort(vs), sorted(np.unique(agg.labels_))))
        order = np.argsort([sortdict[l] for l in agg.labels_])
        
        
        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(mat_raw)}')        
        sns.heatmap(mat_raw[order, :], cmap=cm.bwr, vmin=-3, vmax=3, ax=ax)
        fig.savefig(f'plots/pileups_final/diag50_raw_{tf_name}_heatmap.png')

        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(mat_balanced)}')        
        sns.heatmap(mat_balanced[order, :], cmap=cm.bwr, vmin=-3, vmax=3, ax = ax)
        fig.savefig(f'plots/pileups_final/diag50_balanced_{tf_name}_heatmap.png')

    for tf_name in vals_diag_raw_25:
        mat_raw = deepcopy(vals_diag_raw_25[tf_name])
        mat_balanced = deepcopy(vals_diag_balanced_25[tf_name])
    
        mat_raw[np.isnan(mat_raw)] = 0
        mat_balanced[np.isnan(mat_balanced)] = 0
        
        
        agg = AgglomerativeClustering(n_clusters = 6)
        agg = agg.fit((mat_raw))
        order = np.argsort(agg.labels_)

        vs = []
        for i in range(np.max(agg.labels_)+1):
            vs.append(np.mean(mat_raw[:, 80:120][agg.labels_ == i]))
        sortdict = dict(zip(np.argsort(vs), sorted(np.unique(agg.labels_))))
        order = np.argsort([sortdict[l] for l in agg.labels_])
        
        
        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(mat_raw)}')        
        sns.heatmap(mat_raw[order, :], cmap=cm.bwr, vmin=-3, vmax=3, ax=ax)
        fig.savefig(f'plots/pileups_final/diag25_raw_{tf_name}_heatmap.png')

        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(mat_balanced)}')        
        sns.heatmap(mat_balanced[order, :], cmap=cm.bwr, vmin=-3, vmax=3, ax = ax)
        fig.savefig(f'plots/pileups_final/diag25_balanced_{tf_name}_heatmap.png')
        

    for tf_name in vals_raw:
        mat_raw = deepcopy(vals_raw[tf_name])
        mat_balanced = deepcopy(vals_balanced[tf_name])
        print(tf_name)
        print(mat_raw.shape)
        print(type(mat_raw), type(mat_balanced))

        mat_raw[np.isnan(mat_raw)] = 0
        mat_balanced[np.isnan(mat_balanced)] = 0
        
        agg = AgglomerativeClustering(n_clusters = 6)
        agg = agg.fit((mat_raw))

        order = np.argsort(agg.labels_)
        vs = []
        for i in range(np.max(agg.labels_)+1):
            vs.append(np.mean(mat_raw[:, 72:128][agg.labels_ == i]))
        sortdict = dict(zip(np.argsort(vs), sorted(np.unique(agg.labels_))))
        order = np.argsort([sortdict[l] for l in agg.labels_])
        
        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(mat_raw)}')        
        sns.heatmap(mat_raw[order, :], cmap=cm.bwr, vmin=-3, vmax=3, ax=ax)
        fig.tight_layout()
        fig.savefig(f'plots/pileups_final/crossdiag_raw_{tf_name}_heatmap.png')

        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(mat_balanced)}')
        sns.heatmap(mat_balanced[order, :], cmap=cm.bwr, vmin=-3, vmax=3, ax = ax)
        fig.tight_layout()
        fig.savefig(f'plots/pileups_final/crossdiag_balanced_{tf_name}_heatmap.png')    


        
    for tf_name in mats_raw:       
        pileup_raw = np.nanmean(mats_raw[tf_name], axis=0)
        pileup_balanced = np.nanmean(mats_balanced[tf_name], axis=0)    
        
        pileup_raw[np.isnan(pileup_raw)] = 0
        pileup_balanced[np.isnan(pileup_balanced)] = 0
        
        
        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(pileup_raw)}')
        sns.heatmap(pileup_raw, cmap=cm.bwr, vmin=-1, vmax=1, ax=ax)
        fig.savefig(f'plots/pileups_final/diag_raw_{tf_name}_pileup.png')

        fig, ax = plt.subplots()
        ax.set_title(f'{tf_name}: {len(pileup_balanced)}')
        sns.heatmap(pileup_balanced, cmap=cm.bwr, vmin=-1, vmax=1, ax=ax)    
        fig.savefig(f'plots/pileups_final/diag_balanced_{tf_name}_pileup.png')
        
        
    ctcf_raw = deepcopy(vals_diag_raw_50['ctcf_peaks'])
    foxp3_raw = deepcopy(vals_diag_raw_50['foxp3_peaks_yuri'])
    fig, ax = plt.subplots()
    ax.plot(np.nanmean(ctcf_raw, axis=0), label=f'CTCF: {len(ctcf_raw)}')
    ax.plot(np.nanmean(foxp3_raw, axis=0), label=f'FoxP3: {len(foxp3_raw)}')

    ax.legend()
    fig.savefig(f'plots/pileups_final/plot_diag50.png')


    ctcf_raw = deepcopy(vals_diag_raw_25['ctcf_peaks'])
    foxp3_raw = deepcopy(vals_diag_raw_25['foxp3_peaks_yuri'])
    fig, ax = plt.subplots()
    ax.plot(np.nanmean(ctcf_raw, axis=0), label=f'CTCF: {len(ctcf_raw)}')
    ax.plot(np.nanmean(foxp3_raw, axis=0), label=f'FoxP3: {len(foxp3_raw)}')
    ax.legend()
    fig.savefig(f'plots/pileups_final/plot_diag25.png')


    ctcf_raw = deepcopy(vals_raw['ctcf_peaks'])
    foxp3_raw = deepcopy(vals_raw['foxp3_peaks_yuri'])
    fig, ax = plt.subplots()
    ax.plot(np.nanmean(ctcf_raw, axis=0), label=f'CTCF: {len(ctcf_raw)}')
    ax.plot(np.nanmean(foxp3_raw, axis=0), label=f'FoxP3: {len(foxp3_raw)}')
    ax.legend()
    fig.savefig(f'plots/pileups_final/plot_crossdiag.png')
    return mats_raw

import itertools
def cdf(fcs, names, xlabel="", folder="tmp", filename="tmp", xmin=None, xmax=None):
    
    lens ={}
    if xmin is None:
        xmin = np.min(fcs[0])
    if xmax is None:
        xmax = np.max(fcs[0])
    xs = np.linspace(xmin, xmax, 1000)
    fc_dict = {}
    for i, fc in enumerate(fcs):
        name = names[i]
        fc_dict[name] = []
        for j, x in enumerate(xs):
            fc_dict[name].append((np.asarray(fc) < x).mean())
        lens[name] = len(fc)
    fig, axs = plt.subplots(figsize=(5, 3))
    for name in fc_dict:
        axs.plot(xs, fc_dict[name], label=f'{name}, {lens[name]}')
    axs.legend()

    axs.set_xlabel(xlabel)
    axs.set_xlim([xmin, xmax])
    plt.tight_layout()

    dict_for_ks = dict(zip(names, fcs))
    for pair in itertools.combinations(list(dict_for_ks.keys()), 2):
        name1, name2 = pair
        ks_test = scipy.stats.ks_2samp(dict_for_ks[name1], dict_for_ks[name2])
        pval = ks_test[1] 
        print(pval)

    return fig, axs
    # pval = ks_test[1]
    # axs.set_title(pval.round(3))
    # axs[1].plot(xs, ys_treatment-ys_control, label=f'{name1} - {name2}')





import scipy
import seaborn as sns
def kiko_v_tcon_loop_binding(loop_dict, atac_dict, partners):
    for loop_name in loop_dict:
        loops = loop_dict[loop_name]
        enrichment_df = pd.DataFrame()
        pval_df = pd.DataFrame()

        for atac_name in atac_dict:
            for partner1_name in list(partners.keys()):
                atac_kiko, atac_tcon = atac_dict[atac_name]
                enrichments = []
                pvals = []
                for partner1_name in list(partners.keys()):
                    partner1 = partners[partner1_name]
                    ancs_kiko = (atac_kiko).intersect(loops.slop(b=5000, g='./annotations/chromsizes'), u=True)
                    n_kiko = len(ancs_kiko)
                    success_kiko = ancs_kiko.intersect(partner1, u=True)
                    n_success_kiko = len(success_kiko)
                    probability_kiko = n_success_kiko / n_kiko

                    ancs_tcon = (atac_tcon).intersect(loops.slop(b=5000, g='./annotations/chromsizes'), u=True)
                    n_tcon = len(ancs_tcon)
                    success_tcon = ancs_tcon.intersect(partner1, u=True)
                    n_success_tcon = len(success_tcon)
                    probability_tcon = n_success_tcon / n_tcon


                    if probability_tcon != 0:
                        enrichment = probability_kiko / probability_tcon
                        phat = (n_success_kiko+n_success_tcon)/(n_kiko + n_tcon)
                        zscore = (probability_tcon-probability_kiko)/np.sqrt(phat*(1-phat)*(1/n_kiko + 1/n_tcon))
                        pval = 2*(1-scipy.stats.norm.cdf(np.abs(zscore), 0, 1))    
                    else:
                        enrichment = 0
                        pval = 1
                    enrichments.append(enrichment)
                    pvals.append(pval)
                enrichment_df[atac_name] = enrichments
                enrichment_df.index = list(partners.keys())

                pval_df[atac_name] = pvals
                pval_df.index = list(partners.keys())

        fig, ax = plt.subplots()
        sns.heatmap(np.log2(enrichment_df), annot=True, fmt='.2f',
           vmin=-1, vmax = 1, cmap=cm.bwr, ax=ax)
        fig.tight_layout()
        fig.savefig(f'./plots/kiko_v_tcon_loop_binding/{loop_name}_enrichment.png')


        fig, ax = plt.subplots()
        sns.heatmap((pval_df), annot=True, fmt='.2f',
           vmin=0, vmax = .1, cmap=cm.bwr, ax=ax)
        fig.tight_layout()
        fig.savefig(f'./plots/kiko_v_tcon_loop_binding/{loop_name}_pvalue.png')


def differential_hic_scatterplot(mustache_peaks_loops, mustache_nonsig_loops, cool_treg, cool_tconv):
    tconv_scores = []
    treg_scores = []
    print("Peaks", len(mustache_peaks_loops), "Nonsig", len(mustache_nonsig_loops))
    wsz = 4
    for i in mustache_peaks_loops:
        l1, l2 = i[:3], i[3:6]
        
        l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)
        
        val_tconv = cool_tconv.matrix(balance=True).fetch(l1, l2).sum()
        val_treg = cool_treg.matrix(balance=True).fetch(l1, l2).sum()
        
        tconv_scores.append(val_tconv)
        treg_scores.append(val_treg)
        
    print("Halfway")
    tconv_shifted_scores = []
    treg_shifted_scores = []
    shift = 20
    for i in mustache_nonsig_loops:
        l1, l2 = i[:3], i[3:6]
        l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)
        val_tconv = cool_tconv.matrix(balance=True).fetch(l1, l2).sum()
        val_treg = cool_treg.matrix(balance=True).fetch(l1, l2).sum()
        
        tconv_shifted_scores.append(val_tconv)
        treg_shifted_scores.append(val_treg)
    
    tconv_scores, treg_scores, tconv_shifted_scores, treg_shifted_scores = np.asarray(tconv_scores), np.asarray(treg_scores), np.asarray(tconv_shifted_scores), np.asarray(treg_shifted_scores)


    fig, ax = plt.subplots()
    ax.scatter(np.log2(tconv_scores), np.log2(treg_scores), s=5)
    # ax.set_xlim([3, 13])
    # ax.set_ylim([3, 13])
    ax.set_title("Differential loops")
    ax.set_xlabel("T-reg Hi-C signal (log2)")
    ax.set_ylabel("T-conv Hi-C signal (log2)")
    fig.savefig('./plots/qc/differential_loops.png')


    fig, ax = plt.subplots()
    ax.scatter(np.log2(tconv_shifted_scores), np.log2(treg_shifted_scores), s=5)
    # ax.set_xlim([3, 13])
    # ax.set_ylim([3, 13])
    ax.set_title("Non-differential loops")
    ax.set_xlabel("T-reg Hi-C signal (log2)")
    ax.set_ylabel("T-conv Hi-C signal (log2)")
    fig.savefig('./plots/qc/nondifferential_loops.png')


    fig, ax = plt.subplots()
    ax.scatter(np.log2(tconv_scores), np.log2(treg_scores), c='red', alpha=.5, s=5, label='Real loops (differential)')
    ax.scatter(np.log2(tconv_shifted_scores), np.log2(treg_shifted_scores), c='blue', alpha=.5, s=5, label='Nonsig loops')
    # ax.set_xlim([0, .1])
    # ax.set_ylim([0, .1])
    ax.set_title("Comparison")
    ax.set_xlabel("T-reg Hi-C signal (log2)")
    ax.set_ylabel("T-conv Hi-C signal (log2)")
    fig.savefig('./plots/qc/comparison_loops.png')

    plt.legend()     

def binding_site_anchor_size(full_anchors, tfs, res=5000, shift=50, suffix=""):
    frac_dict = {}
    shifted_frac_dict = {}
    for tf_name in tfs:
        fracs = []
        shifted_fracs = []
        tf = tfs[tf_name]
        xs = []
        for i in [0, 1, 2, 3, 4, 5, 8, 10, 15, 20]:
            i = i*res
            frac = len(tf.intersect(full_anchors.slop(b=i, g='./annotations/chromsizes'), u=True))/len(tf)
            shifted_frac = len(tf.intersect(full_anchors.shift(s=shift*res, g='./annotations/chromsizes').slop(b=i, g='./annotations/chromsizes'), u=True))/len(tf)    
            xs.append(i)
            fracs.append(frac)
            shifted_fracs.append(shifted_frac)
        frac_dict[tf_name]= fracs
        shifted_frac_dict[tf_name] = shifted_fracs
    xs = np.asarray(xs)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))   
    for tf_name in frac_dict:
        axs[0].plot(xs/5000, frac_dict[tf_name], label=tf_name, marker='o')
    axs[0].legend()   
    for tf_name in frac_dict:
        axs[1].plot(xs/5000, shifted_frac_dict[tf_name], label=tf_name, marker='o')
    axs[1].legend()
    axs[0].set_ylim([0, 1])
    axs[1].set_ylim([0, 1])
    fig.savefig(f'./plots/qc/fraction_TF_covered_{suffix}.png')


    frac_dict = {}
    shifted_frac_dict = {}
    for tf_name in tfs:
        fracs = []
        shifted_fracs = []
        tf = tfs[tf_name]
        xs = []
        for i in [0, 1, 2, 3, 4, 5, 8, 10, 15, 20]:
            i = i*res
            frac = len(full_anchors.slop(b=i, g='./annotations/chromsizes').intersect(tf, u=True))/len(full_anchors)
            shifted_frac = len(full_anchors.shift(s=shift*res, g='./annotations/chromsizes').slop(b=i, g='./annotations/chromsizes').intersect(tf, u=True))/len(full_anchors)    
            xs.append(i)
            fracs.append(frac)
            shifted_fracs.append(shifted_frac)
        frac_dict[tf_name]= fracs
        shifted_frac_dict[tf_name] = shifted_fracs
    xs = np.asarray(xs)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))   
    for tf_name in frac_dict:
        axs[0].plot(xs/5000, frac_dict[tf_name], label=tf_name, marker='o')
    axs[0].legend()   
    for tf_name in frac_dict:
        axs[1].plot(xs/5000, shifted_frac_dict[tf_name], label=tf_name, marker='o')
    axs[1].legend()
    axs[0].set_ylim([0, 1])
    axs[1].set_ylim([0, 1])    
    fig.savefig(f'./plots/qc/fraction_LOOPS_covered_{suffix}.png')    

import time
def apa(loops, cool_treg, cool_tconv, expected_treg, expected_tconv, wsz=20, name='loops', balance=False, pc = .1):
    mats_treg = []
    df = expected_treg
    df['region'] = df['region'].apply(lambda x: str(x))
    c = 0
    for i in sorted(list(loops)):
        t2 = time.time()

        l1, l2 = i[:3], i[3:6]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l2[1])
        if (e-s)/5000 < 2*wsz:
            continue
        t1 = time.time()
        expected = make_expected(df, chrom, s, e, wsz+1, res=5000, balance=balance)
        t1f = time.time()
        #print(t1f-t1)
        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        new_l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)

        val = np.asarray(cool_treg.matrix(balance=balance).fetch(new_l1, new_l2))

        if val.shape != (2*wsz+3, 2*wsz+3):
            print(l1, l2)
            print(val.shape)
            continue
        mats_treg.append(np.log2((val+expected*pc)/(expected+expected*pc)))
        c += 1
        #t2f = time.time()

    mats_tconv = []
    df = expected_tconv
    df['region'] = df['region'].apply(lambda x: str(x))
    places = []
    for i in sorted(list(loops)):
        l1, l2 = i[:3], i[3:6]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l2[1])
        if (e-s)/5000 < 2*wsz:
            #print((e-s)/5000)
            continue
        expected = make_expected(df, chrom, s, e, wsz+1, res=5000, balance=balance)

        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        new_l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)

        val = np.asarray(cool_tconv.matrix(balance=balance).fetch(new_l1, new_l2))
        if val.shape != (2*wsz+3, 2*wsz+3):
            print(l1, l2)
            continue
        places.append(l1 + l2)
        mats_tconv.append(np.log2((val+expected*pc)/(expected+expected*pc)))
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    c = axs[0].matshow(np.nanmean(mats_treg, axis=0), vmin=-3, vmax=3, cmap=cm.bwr)
    axs[1].matshow(np.nanmean(mats_tconv, axis=0), vmin=-3, vmax=3, cmap=cm.bwr)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(c, cax=cbar_ax)
    fig.savefig(f'./plots/qc/{name}.png')
    return np.asarray(mats_treg), np.asarray(mats_tconv), places

def apa_remove_zero(loops, cool_treg, cool_tconv, expected_treg, expected_tconv, wsz=20, name='loops', balance=False, pc = .1):
    mats_treg = []
    df = expected_treg
    df['region'] = df['region'].apply(lambda x: str(x))
    c = 0
    for i in loops:
        t2 = time.time()

        l1, l2 = i[:3], i[3:6]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l2[1])
        if (e-s)/5000 < 2*wsz:
            continue
        t1 = time.time()
        expected = make_expected(df, chrom, s, e, wsz+1, res=5000, balance=balance)
        t1f = time.time()
        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        new_l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)

        val = np.asarray(cool_treg.matrix(balance=balance).fetch(new_l1, new_l2))
        val[val==0] = np.nan
        if val.shape != (2*wsz+3, 2*wsz+3):
            print(l1, l2)
            print(val.shape)
            continue
        mats_treg.append(np.log2((val+expected*pc)/(expected+expected*pc)))
        c += 1

    mats_tconv = []
    df = expected_tconv
    df['region'] = df['region'].apply(lambda x: str(x))
    for i in loops:
        break
        l1, l2 = i[:3], i[3:6]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l2[1])
        if (e-s)/5000 < 2*wsz:
            continue
        expected = make_expected(df, chrom, s, e, wsz+1, res=5000, balance=balance)

        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        new_l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)

        val = np.asarray(cool_tconv.matrix(balance=balance).fetch(new_l1, new_l2))
        val[val==0] = np.nan
        if val.shape != (2*sz+3, 2*wsz+3):
            print(l1, l2)
            continue
        mats_tconv.append(np.log2((val+expected*pc)/(expected+expected*pc)))
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    c = axs[0].matshow(np.nanmean(mats_treg, axis=0), vmin=-3, vmax=3, cmap=cm.bwr)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(c, cax=cbar_ax)
    fig.savefig(f'./plots/qc/{name}.png')
    return np.asarray(mats_treg), np.asarray(mats_tconv)


def apa_diag(anchors, cool_treg, cool_tconv, expected_treg, expected_tconv, balance=True, wsz=20, name='loops', pc_ratio=.1):
    mats_treg = []
    df = expected_treg
    df['region'] = df['region'].apply(lambda x: str(x))
    c = 0
    for i in anchors:
        l1 = i[:3]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l1[1])

        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)

        try:
            val = np.asarray(cool_treg.matrix(balance=balance).fetch(new_l1))
        except Exception as e:
            if "Genomic region out of bounds" in str(e):
                continue
            else:
                break

        expected = make_diag_expected(df, chrom, 2*wsz+3, res=5000, balance=balance)

        if val.shape != (2*wsz+3, 2*wsz+3):
            print(l1, l1)
            print(val.shape)
            continue

        treg_ratio = np.log2((val+expected*pc_ratio)/(expected+expected*pc_ratio))
        treg_ratio[np.isinf(treg_ratio)] = 0
        mats_treg.append(treg_ratio)
        c += 1

    mats_tconv = []
    df = expected_tconv
    df['region'] = df['region'].apply(lambda x: str(x))
    for i in anchors:
        l1 = i[:3]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l1[1])


        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)

        try:
            val = np.asarray(cool_tconv.matrix(balance=balance).fetch(new_l1))
        except Exception as e:
            if "Genomic region out of bounds" in str(e):
                continue
            else:
                break

        expected = make_diag_expected(df, chrom, 2*wsz+3, res=5000, balance=balance)

        if val.shape != (2*wsz+3, 2*wsz+3):
            print(l1, l2)
            continue

        tcon_ratio = (np.log2((val+expected*pc_ratio)/(expected+expected*pc_ratio)))
        tcon_ratio[np.isinf(tcon_ratio)] = 0
        mats_tconv.append(tcon_ratio)

    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    c = axs[0].matshow(np.nanmean(mats_treg, axis=0), vmin=-1, vmax=1, cmap=cm.bwr)
    axs[1].matshow(np.nanmean(mats_tconv, axis=0), vmin=-1, vmax=1, cmap=cm.bwr)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(c, cax=cbar_ax)
    fig.savefig(f'./plots/qc/{name}.png')
    return np.asarray(mats_treg), np.asarray(mats_tconv)





def make_diag_expected(df, chrom, wsz, res=5000, balance=False):
    chrom = str(chrom)
    expected = np.zeros((wsz, wsz))
    n = expected.shape[0]

    test = np.zeros((wsz, wsz))
    n = test.shape[0]
    averages = df[(df['region'] == chrom)].sort_values('diag')
    for i, place in enumerate(range(2, wsz)):
        if balance==True:
            expected_val = averages.iloc[place, :]['balanced.avg']
        else:
            expected_val = averages.iloc[place, :]['count.avg']
        np.fill_diagonal(expected[place:, :], expected_val)
    result = expected+expected.T - np.diag(np.diag(expected))
    result[result==0]=np.nan
    return result



def apa_hicdc(loops, cool_treg, cool_tconv, expected_treg, expected_tconv, wsz=20, name='loops'):
    mats_treg = []
    df = expected_treg
    df['region'] = df['region'].apply(lambda x: str(x))    
    c = 0
    for i in loops:
        l1, l2 = i[:3], i[3:6]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l2[1])
        if (e-s)/5000 < 2*wsz:
            continue
        expected = make_expected(df, chrom, s, e, wsz, res=5000)
        
        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        new_l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)

        val = np.asarray(cool_treg.matrix(balance=False).fetch(new_l1, new_l2))
        if val.shape != (2*wsz+1, 2*wsz+1):
            print(l1, l2)
            pritn(val.shape)
            continue
        mats_treg.append(np.log2((val+1)/(expected+1)))
        c += 1

    mats_tconv = []
    df = expected_tconv
    df['region'] = df['region'].apply(lambda x: str(x))        
    for i in loops:
        l1, l2 = i[:3], i[3:6]
        chrom = l1[0]
        s = int(l1[1])
        e = int(l2[1])
        if (e-s)/5000 < 2*wsz:
            continue
        expected = make_expected(df, chrom, s, e, wsz, res=5000)
        
        new_l1 = (l1[0], int(l1[1])-5000*wsz, int(l1[2])+5000*wsz)
        new_l2 = (l2[0], int(l2[1])-5000*wsz, int(l2[2])+5000*wsz)

        val = np.asarray(cool_tconv.matrix(balance=False).fetch(new_l1, new_l2))
        if val.shape != (2*wsz+1, 2*wsz+1):
            print(l1, l2)
            continue
        mats_tconv.append(np.log2((val+1)/(expected+1)))
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    c = axs[0].matshow(np.nanmean(mats_treg, axis=0), vmin=-3, vmax=3, cmap=cm.bwr)
    axs[1].matshow(np.nanmean(mats_tconv, axis=0), vmin=-3, vmax=3, cmap=cm.bwr)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(c, cax=cbar_ax)
    fig.savefig(f'./plots/qc/hicdc/{name}.png')    
    return np.asarray(mats_treg), np.asarray(mats_tconv)


def shifted_tf_binding(partners, mustache_tcon_anchors, mustache_treg_anchors, mustache_nonsig_anchors):
    for name in partners:
        poop = partners[name]

        xs, ys = [], []
        xs_2, ys_2 = [], []
        xs_3, ys_3 = [], []
        for i in range(-20, 20, 2):
            shift = i
            places = mustache_tcon_anchors.shift(s=i*5000, g='./annotations/chromsizes').intersect(poop, u=True)
            x = shift
            y = len(places)
            xs.append(x)
            ys.append(y)

            places2 = mustache_treg_anchors.shift(s=i*5000, g='./annotations/chromsizes').intersect(poop, u=True)
            xs_2.append(shift)
            ys_2.append(len(places2))

            places3 = mustache_nonsig_anchors.shift(s=i*5000, g='./annotations/chromsizes').intersect(poop, u=True)
            xs_3.append(shift)
            ys_3.append(len(places3))
            
        [xs, ys, xs_2, ys_2, xs_3,  ys_3] = map(np.asarray, [xs, ys, xs_2, ys_2, xs_3, ys_3])
        fig, axs = plt.subplots( figsize=(6, 4))
        axs.set_xlabel('Scale = 5kb')
        axs.plot(xs, ys/len(mustache_tcon_anchors), marker='o', label='Tcon anchors')
        axs.plot(xs_2, ys_2/len(mustache_treg_anchors), marker='o', label='Treg anchors')
        axs.plot(xs_3, ys_3/len(mustache_nonsig_anchors), marker='o', label='Non-sig anchors')
        axs.set_title(name)
        # ax.set_title('Anchors')
        axs.legend()
        fig.savefig(f'plots/shifted_tf_binding/{name}.png')


import time 
def make_expected(df, chrom, s, e, wsz, res=5000, balance=False):
    chrom = str(chrom)
    expected = np.zeros((wsz*2+1, wsz*2+1))
    n = expected.shape[0]

    test = np.zeros((wsz*2+1, wsz*2+1))
    n = test.shape[0]
    diag = (e-s)//res
    averages = df[(df['region'] == chrom)].sort_values('diag')
    for i, place in enumerate(range(-2*wsz-1, 2*wsz+1)):
        if balance==True:
            expected_val = averages.iloc[diag+place, :]['balanced.avg']
        else:
            expected_val = averages.iloc[diag+place, :]['count.avg']
        if place < 0:
            np.fill_diagonal(expected[-place:, :], expected_val)        
        if place >=0:
            np.fill_diagonal(expected[:, place:], expected_val)        

    return expected


import seaborn as sns
import itertools
from copy import deepcopy

def atac_vs_expected(anclist, names, partners, atacs):
        # enrichment_df.T.plot(ax=ax, linewidth=0, marker='o')
        # fig.savefig(f'./plots/atac_vs_expected_heatmap/{names[i]}.png')
        # ax.set_ylim(-1.5, 1.5)
    xs, ys, hues = [], [], []

    hypergeom_ps = []
    for i, ancs_of_interest in enumerate(anclist):
        pvals = []
        pval_loc = []
        enrichment_df = pd.DataFrame()
        for partner1_name in list(partners.keys()):
            enrichments = []
            for atac_name in list(atacs.keys()):

                tf = partners[partner1_name]
                atac = atacs[atac_name]
                cross_binders = 0
                default = atac.intersect(tf, u=True)
                n_success_genomwide = len(default)
                n_genomewide = len(atac)
                genomwide_p = len(default)/len(atac)

                atac_in_ancs = atac.intersect(ancs_of_interest, u=True)
                atac_in_ancs_with_tf = atac_in_ancs.intersect(tf, u=True)
                anchor_p = len(atac_in_ancs_with_tf) / len(atac_in_ancs)
                n_success_ancs = len(atac_in_ancs_with_tf)
                n_in_ancs = len(atac_in_ancs)
                if genomwide_p != 0:
                    enrichment = np.log2(anchor_p / genomwide_p)
                    phat = (n_success_ancs+n_success_genomwide)/(n_genomewide + n_in_ancs)
                    zscore = (anchor_p-genomwide_p)/np.sqrt(phat*(1-phat)*(1/n_genomewide + 1/n_in_ancs))
                    pval = (1-scipy.stats.norm.cdf(np.abs(zscore), 0, 1))  
                    pvals.append(pval)
                    pval_loc.append(partner1_name)
                else:
                    enrichment = 1e-5
                if i == 0:
                    xs.append(partner1_name)
                    ys.append(genomwide_p)
                    hues.append('Genome-wide')


                xs.append(partner1_name)
                ys.append(anchor_p)
                hues.append(names[i])
                enrichments.append(enrichment)

            enrichment_df[partner1_name] = enrichments

    poop = pd.DataFrame()
    poop['TF'] = xs
    poop['% ATAC overlapping TF'] = ys
    poop['Type'] = hues

    
    M = len(atacs['All ATAC peaks'])
    for tf in np.unique(xs):
        for i in range(0, len(names)):
            anctype = names[i]
            K = len(atacs['All ATAC peaks'].intersect(anclist[i], u=True))
            prop_genomewide = poop[(poop['TF']==tf) & (poop['Type'] == 'Genome-wide')]['% ATAC overlapping TF'].values[0]
            prop_anc = poop[(poop['TF']==tf) & (poop['Type'] == anctype)]['% ATAC overlapping TF'].values[0]
            m = int(prop_genomewide*M)
            k = int(prop_anc*K)
            # print("___", poop[(poop['TF']==tf) & (poop['Type'] == 'Genome-wide')]['% ATAC overlapping TF'])
            print(anctype, tf, prop_genomewide, prop_anc)
            print(scipy.stats.hypergeom.cdf(k, M, m, K))
            print(k, M, m, K)

    fig = sns.catplot(x='TF', y='% ATAC overlapping TF', hue='Type', data=poop, kind='bar', height=4.27, aspect=11.7/7.27)
    fig.set(ylim=(0, .14))

    # for c, _ in enumerate(pvals):
    #     print(names[i], list(partners.keys())[c], pvals[c])
    fig.savefig(f'./plots/atac_vs_expected_heatmap/{names[i]}')


def binding_bars(anclist, names, partners, atacs):
        # enrichment_df.T.plot(ax=ax, linewidth=0, marker='o')
        # fig.savefig(f'./plots/atac_vs_expected_heatmap/{names[i]}.png')
        # ax.set_ylim(-1.5, 1.5)
    xs, ys, hues = [], [], []

    hypergeom_ps = []
    for i, ancs_of_interest in enumerate(anclist):
        pvals = []
        pval_loc = []
        enrichment_df = pd.DataFrame()
        for partner1_name in list(partners.keys()):
            enrichments = []
            for atac_name in list(atacs.keys()):

                tf = partners[partner1_name]
                atac = atacs[atac_name]
                cross_binders = 0
                default = atac.intersect(tf, u=True)
                n_success_genomwide = len(default)
                n_genomewide = len(atac)
                genomwide_p = len(default)/len(atac)

                atac_in_ancs = atac.intersect(ancs_of_interest, u=True)
                atac_in_ancs_with_tf = atac_in_ancs.intersect(tf, u=True)
                anchor_p = len(atac_in_ancs_with_tf) / len(atac_in_ancs)
                n_success_ancs = len(atac_in_ancs_with_tf)
                n_in_ancs = len(atac_in_ancs)
                if genomwide_p != 0:
                    enrichment = np.log2(anchor_p / genomwide_p)
                    phat = (n_success_ancs+n_success_genomwide)/(n_genomewide + n_in_ancs)
                    zscore = (anchor_p-genomwide_p)/np.sqrt(phat*(1-phat)*(1/n_genomewide + 1/n_in_ancs))
                    pval = (1-scipy.stats.norm.cdf(np.abs(zscore), 0, 1))  
                    pvals.append(pval)
                    pval_loc.append(partner1_name)
                else:
                    enrichment = 1e-5
                if i == 0:
                    xs.append(partner1_name)
                    ys.append(genomwide_p)
                    hues.append('Genome-wide')


                xs.append(partner1_name)
                ys.append(anchor_p)
                hues.append(names[i])
                enrichments.append(enrichment)

            enrichment_df[partner1_name] = enrichments

    poop = pd.DataFrame()
    poop['TF'] = xs
    poop['Fraction of ATAC peaks overlapping...'] = ys
    poop['Type'] = hues

    
    M = len(atacs['All ATAC peaks'])
    for tf in np.unique(xs):
        for i in range(0, len(names)):
            anctype = names[i]
            K = len(atacs['All ATAC peaks'].intersect(anclist[i], u=True))
            prop_genomewide = poop[(poop['TF']==tf) & (poop['Type'] == 'Genome-wide')]['Fraction of ATAC peaks overlapping...'].values[0]
            prop_anc = poop[(poop['TF']==tf) & (poop['Type'] == anctype)]['Fraction of ATAC peaks overlapping...'].values[0]
            m = int(prop_genomewide*M)
            k = int(prop_anc*K)
            # print("___", poop[(poop['TF']==tf) & (poop['Type'] == 'Genome-wide')]['% ATAC overlapping TF'])
            print(anctype, tf, prop_genomewide, prop_anc, scipy.stats.hypergeom.cdf(k, M, m, K))

    fig = sns.catplot(x='TF', y='Fraction of ATAC peaks overlapping...', hue='Type', data=poop, kind='bar', height=4.27, aspect=11.7/7.27)

    # for c, _ in enumerate(pvals):
    #     print(names[i], list(partners.keys())[c], pvals[c])
    fig.savefig(f'./plots/binding_bars/{names[i]}')


def foxp3_ctcf_distance(mustache_full_loops, ctcf_peaks, foxp3_peaks,):
    dists_ctcf = []
    seen = {}
    for i in mustache_full_loops.pair_to_bed(ctcf_peaks, type='both'):
        place = tuple(i[:6])
        if seen.get(place):
            continue
        else:
            seen[place]=1
        s1, s2, = i[2], i[5]
        s1, s2 = map(int, [s1, s2])
        dists_ctcf.append((s2-s1)/5000)

        
    dists_foxp3 = []
    seen = {}
    for i in mustache_full_loops.pair_to_bed(foxp3_peaks, type='both'):
        place = tuple(i[:6])
        if seen.get(place):
            continue
        else:
            seen[place]=1
        s1, s2, = i[2], i[5]
        s1, s2 = map(int, [s1, s2])
        dists_foxp3.append((s2-s1)/5000)

    conditions = {
    'foxp3':dists_foxp3,
    'ctcf':dists_ctcf,
    }
    cdf([dists_foxp3, dists_ctcf], ['FoxP3', 'CTCF'], xlabel='Bins separating loop (5kb)',
       folder='distance_cdf', filename='foxp3_vs_ctcf.png') 

import random
def make_df_nice(df):
    new_df = deepcopy(df)
    new_df[np.isnan(new_df)] = 0
    new_df[np.isinf(new_df)] = 0
    new_df.iloc[:, 0] += 1
    return new_df

def make_order2_custom(matrix, method='average', metric='cosine'):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method=method, metric=metric)
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)

        order = dendro['leaves']
    else:
        order = np.arange(len(matrix))
    # order = dendro['leaves']
    return order        


def make_order2(matrix):
    if len(matrix) > 1:
        linkage = scipy.cluster.hierarchy.linkage(matrix, method='average', metric='cosine')
        dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                            color_threshold=-np.inf)

        order = dendro['leaves']
    else:
        order = np.arange(len(matrix))
    # order = dendro['leaves']
    return order        

def both_looping_heatmap(looplist, names, partners, ancsets, prefix=''):
    df_dict = {}
    for i, (name, loops_of_interest) in enumerate(looplist.items()):
        listified = list(loops_of_interest)
        left_anchors = []
        right_anchors = []
        all_anchors = []
        for x in listified:
            l1, l2 = x[:3], x[3:6]
            left_anchors.append(tuple(l1))
            right_anchors.append(tuple(l2))
        all_anchors = left_anchors + right_anchors
        left_anchors = np.asarray(left_anchors)
        right_anchors = np.asarray(right_anchors)
        fc_enrichment_df = pd.DataFrame()
        subtraction_enrichment_df = pd.DataFrame()
        pvalue_df = pd.DataFrame()
        pvalue2_df = pd.DataFrame()

        for partner1_name in list(partners.keys()):
            fc_enrichments = []
            subtraction_enrichments = []
            pvals = []
            pvals2 = []
            for partner2_name in list(partners.keys()):
                if partner1_name < partner2_name:
                    fc_enrichments.append(0)
                    subtraction_enrichments.append(0)
                    pvals.append(0)
                    pvals2.append(0)
                    continue

                partner1 = partners[partner1_name]
                partner2 = partners[partner2_name]
                cross_binders = 0
                for loop in listified:
                    l1 = tuple(loop[:3])
                    l2 = tuple(loop[3:6])
                    if ((l1 in ancsets[partner1_name]) and (l2 in ancsets[partner2_name])) or \
                        ((l2 in ancsets[partner1_name]) and (l1 in ancsets[partner2_name])):
                            cross_binders += 1
                t1 = time.time()
                n = 100
                binding_prop = []
                for k in range(n):
                    permutation = np.random.permutation(len(left_anchors))
                    left_anchors = left_anchors[permutation]
                    rand_cross_binders = 0
                    for _, l1 in enumerate(right_anchors):
                        l2 = left_anchors[_]
                        if ((tuple(l1) in ancsets[partner1_name]) and (tuple(l2) in ancsets[partner2_name])) or \
                            ((tuple(l2) in ancsets[partner1_name]) and (tuple(l1) in ancsets[partner2_name])):
                                rand_cross_binders += 1
                    binding_prop.append(rand_cross_binders/len(listified))
                loop_p = cross_binders/len(listified)
                rand_ps = np.asarray(binding_prop)

                print("Effect", np.min(rand_ps), np.max(rand_ps), loop_p, cross_binders/len(listified), partner1_name, partner2_name)
                # pval2 = min((rand_ps < loop_p).mean(), (rand_ps > loop_p).mean())
                # pvals2.append(pval2)

                # n = 2_000
                # rand_cross_binders = 0
                # for _ in range(n):
                #     randloop = random.sample(all_anchors, 2)
                #     l1, l2 = randloop
                #     l1, l2 = tuple(l1[:3]), tuple(l2[:3])
                #     if ((l1 in ancsets[partner1_name]) and (l2 in ancsets[partner2_name])) or \
                #         ((l2 in ancsets[partner1_name]) and (l1 in ancsets[partner2_name])):
                #             rand_cross_binders += 1
                # print(rand_cross_binders/n)
                # pval = one_sample_proportion_significance_test(cross_binders, len(listified), rand_cross_binders/n)
                # pval = 1

                rand_bothp = np.mean(binding_prop)
                real_bothp = cross_binders/len(listified)
                enrichment = real_bothp/ rand_bothp
                subtraction_enrichment = real_bothp - rand_bothp
                fc_enrichments.append(enrichment)  
                subtraction_enrichments.append(subtraction_enrichment)   

                print(len(fc_enrichments))
                print(len(subtraction_enrichments))

            fc_enrichment_df[partner1_name] = fc_enrichments
            subtraction_enrichment_df[partner1_name] = subtraction_enrichments
        fc_enrichment_df.index = list(partners.keys())
        subtraction_enrichment_df.index = list(partners.keys())

        fc_enrichment_df += fc_enrichment_df.T - np.diag(np.diag(fc_enrichment_df))
        subtraction_enrichment_df += subtraction_enrichment_df.T - np.diag(np.diag(subtraction_enrichment_df))

        df_dict[f'{prefix}_{name}'] = fc_enrichment_df
        sns.set(font_scale=1)

        if i == 0:
            order = make_order2(make_df_nice(fc_enrichment_df))

        df = fc_enrichment_df.iloc[order, order]
        fig = sns.clustermap(np.log2(df),
           vmin=-1.5, vmax = 1.5, cmap=cm.bwr, row_cluster=False, col_cluster=False, annot=False, figsize=(15, 15), yticklabels=True, xticklabels=True)
        fig.cax.set_visible(False)
        fig.savefig(f'./plots/both_looping_heatmap/{prefix}_{name}.png')        

        subtraction_enrichment_df = subtraction_enrichment_df.iloc[order, order]
        sub_fig = sns.clustermap(subtraction_enrichment_df,
           vmin=-.2, vmax = .2, cmap=cm.bwr, row_cluster=False, col_cluster=False, annot=False, figsize=(15, 15), yticklabels=True, xticklabels=True)
        sub_fig.cax.set_visible(False)
        plt.tight_layout()
        sub_fig.savefig(f'./plots/both_looping_heatmap/{prefix}_subtraction_{name}.png')        
    return df_dict




def all_atac_figures(loop_dict, mustache_full_anchors, all_tcon_v_treg, all_kiko_v_treg,):
    means = {}
    deseq_fcs_tcon = {}

    sizes = {}
    places = {}
    deseq_fcs_kiko = {}

    ancs = mustache_full_anchors

    for i in ancs.intersect(all_tcon_v_treg, wo=True):
        l1 = tuple(i[:3])
        place = l1
        deseq_fcs_tcon.setdefault(place, [])
        deseq_fcs_tcon[place].append(float(i[7]))
        
    for i in ancs.intersect(all_kiko_v_treg, wo=True):
        l1 = tuple(i[:3])
        place = l1
        deseq_fcs_kiko.setdefault(place, [])    
        deseq_fcs_kiko[place].append(float(i[7]))

    poops = {}
    for loop_name in loop_dict:
        loops_of_interest = loop_dict[loop_name]
        xs_tcon = []
        xs_treg = []
        xs_kiko = []
        ys_tcon = []
        ys_treg = []
        ys_kiko = []

        cs = []

        seen = {}
        messup = 0
        for i in loops_of_interest.pair_to_bed(all_tcon_v_treg, type='either'):

            l1, l2 = tuple(i[:3]), tuple(i[3:6])
            place = (l1, l2)
            if seen.get(place, 0):
                continue
            seen[(l1, l2)] = 1
            c = float(i[11])
            for i, l in enumerate([l1, l2]):


                if deseq_fcs_tcon.get(l):
                    pass
                    tcon_fc = np.mean(np.asarray(deseq_fcs_tcon[l]))
                    kiko_fc = np.mean(np.asarray(deseq_fcs_kiko[l]))
                else:
                    tcon_fc = 0
                    kiko_fc = 0
                    messup+= 1
                if i == 0:
                    xs_tcon.append(tcon_fc)
                    xs_kiko.append(kiko_fc)
                if i == 1:
                    ys_tcon.append(tcon_fc)
                    ys_kiko.append(kiko_fc)


                    cs.append(c)
        xs_tcon, ys_tcon, xs_treg, ys_treg, xs_kiko, ys_kiko, cs = map(np.asarray, [xs_tcon, ys_tcon, xs_treg, ys_treg, xs_kiko, ys_kiko, cs])
        
        poop = atac_figures(xs_tcon, ys_tcon, xs_kiko, ys_kiko, cs, loop_name)
        
        where = np.abs(cs) > 1
        atac_figures(xs_tcon[where], ys_tcon[where], xs_kiko[where], ys_kiko[where], cs[where], f'big_fc_{loop_name}')    
        
        where = np.abs(cs) < 1
        atac_figures(xs_tcon[where], ys_tcon[where], xs_kiko[where], ys_kiko[where], cs[where], f'small_fc_{loop_name}')    
        poops[loop_name] = poop
    atac_seq_bars(poops, 'all.png')


def atac_cdf_figures(mustache_full_anchors, anc_dict, all_tcon_v_treg, all_kiko_v_treg):
    deseq_fcs_tcon = []
    deseq_fcs_kiko = []

    for i in all_tcon_v_treg.intersect(mustache_full_anchors, u=True):
        l1 = tuple(i[:3])
        place = l1
        deseq_fcs_tcon.append(float(i[4]))
        
    for i in all_kiko_v_treg.intersect(mustache_full_anchors, u=True):
        l1 = tuple(i[:3])
        place = l1
        deseq_fcs_kiko.append(float(i[4]))
        
    poops = {}
    names = []
    fcs_tcon = []
    fcs_kiko = []
    for anc_name in anc_dict:
        tcon_fcs_in_ancs = []
        kiko_fcs_in_ancs = []

        ancs_of_interest = anc_dict[anc_name]

        messup = 0
        for i in all_tcon_v_treg.intersect(ancs_of_interest, u=True):
            tcon_fcs_in_ancs.append(float(i[4]))
        
        for i in all_kiko_v_treg.intersect(ancs_of_interest, u=True):
            kiko_fcs_in_ancs.append(float(i[4]))
        # return (deseq_fcs_tcon, deseq_fcs_kiko, tcon_fcs_in_ancs, kiko_fcs_in_ancs)
        names.append(anc_name)
        fcs_tcon.append(tcon_fcs_in_ancs)
        fcs_kiko.append(kiko_fcs_in_ancs)
    
    cdf(fcs_tcon, names, folder='atac_cdf', filename=f'tcon_anchors.png', 
        xlabel = 'LFC, Tcon ÷ Treg', xmin=-2, xmax=2)
    
    cdf(fcs_kiko, names, folder='atac_cdf', filename=f'kiko_anchors.png', 
        xlabel = 'LFC, KIKO ÷ Treg', xmin=-2, xmax=2)









# %matplotlib inline
# pvals = []
# log2fcs = []
# conditions = []
# seen_genes = {}
# places = []

# cs = []
# for i in peak_genes.intersect(cds, u=True).filter(
#         lambda x: (float(x[4]) < .05)):
#     name = i[3]
#     if seen_genes.get(name):
#         continue
#     seen_genes[name] = 1
#     padj, log2fc = float(i[4]), float(i[6])

#     pvals.append(padj)
#     log2fcs.append(log2fc)
#     if log2fc < 0:
#         place = tuple(i[:3])
#         places.append(place)
    
# conditions, pvals, log2fcs = map(np.asarray, [conditions, pvals, log2fcs])
# len(pvals)

# fig, ax = plt.subplots()
# c = ax.scatter(log2fcs, -np.log10(pvals), cmap = cm.bwr, vmax=-2, vmin=2, s=4)
# ax.set_xlabel("Log2(FC): Higher = Tconv")
# ax.set_ylabel("-log10(Pvalue)")
# fig.colorbar(c)
# right, left = (log2fcs > 0), (log2fcs < 0)
# ax.text(0.9, 0.9, f'{right.sum().round(2)} \n {right.mean().round(2)}', size=10, color='purple', transform=ax.transAxes)
# ax.text(0.1, 0.9, f'{left.sum().round(2)} \n {left.mean().round(2)}', size=10, color='purple', transform=ax.transAxes)
# ax.set_xlim([-6, 6])

# fig, ax = plt.subplots()
# c = ax.scatter(log2fcs, pvals, cmap = cm.bwr, vmax=-2, vmin=2, s=4)
# ax.set_xlabel("Log2(FC): Higher = Tconv")
# ax.set_ylabel("Pvalue")
# fig.colorbar(c)

# right, left = (log2fcs > 0), (log2fcs < 0)
# ax.text(0.9, 0.9, f'{right.sum().round(2)} \n {right.mean().round(2)}', size=10, color='purple', transform=ax.transAxes)
# ax.text(0.1, 0.9, f'{left.sum().round(2)} \n {left.mean().round(2)}', size=10, color='purple', transform=ax.transAxes)
# ax.set_xlim([-6, 6])


def gene_expression_figures(genes, default, anc_dict, gene_fcs, gene_sizes):
    default_fcs = []
    default_pvals = []
    seen_genes = {}
    for i in genes.intersect(default, u=True):
        pval = float(i[4])
        l2fc = float(i[6])
        name = i[3]
        if seen_genes.get(name):
            continue
        else:
            seen_genes[name] = 1

        default_fcs.append(l2fc)
        default_pvals.append(pval)

    sig_props = []
    names = []
    fcs = []
    dict_for_pval = {}
    for anc_name in anc_dict:
        pvals = []
        log2fcs = []
        seen_genes = {}
        places = []

        ancs = anc_dict[anc_name]
        for i in genes.intersect(ancs, wo=True):
            name = i[3]
            if anc_name == 'Non-diff anchors':
                pass
            elif anc_name == 'Treg anchors':
                if (gene_fcs[name][np.argmin(gene_sizes[name])] > 0 ):
                    continue
            elif anc_name == 'Tcon anchors':
                if (gene_fcs[name][np.argmin(gene_sizes[name])] < 0 ):
                    continue
            else:
                print(anc_name)
                raise Exception


            padj, log2fc = float(i[4]), float(i[6])

            if seen_genes.get(name):
                continue
            else:
                seen_genes[name] = 1
            pvals.append(padj)
            log2fcs.append(log2fc)
    
        fcs.append(log2fcs)
        sig_prop = (np.asarray(pvals) < .05).mean()
        sig_props.append(sig_prop)
        names.append(anc_name)
        default_sig_prop = (np.asarray(default_pvals) < .05).mean()
        
        dict_for_pval[anc_name] = ((np.asarray(pvals)<.05).sum(), len(pvals))

    
    cdf(fcs, names, folder='gene_cdf', filename=f'without_all_genes.png', xlabel='Gene expression LFC, Tcon ÷ Treg')

    fig, ax = plt.subplots()

    seen_genes = {}
    all_pvals = []
    all_l2fcs = []

    for i in genes:
        pval = float(i[4])
        l2fc = float(i[6])
        name = i[3]
        if seen_genes.get(name):
            continue
        else:
            seen_genes[name] = 1

        all_pvals.append(pval)
        all_l2fcs.append(l2fc)
    all_pvals = np.asarray(all_pvals)
    all_l2fcs = np.asarray(all_l2fcs)

    dict_for_pval['All genes'] = ((all_pvals < .05).sum(), len(all_pvals))

    cdf(fcs + [all_l2fcs], names + ['All genes'], folder='gene_cdf', filename=f'with_all_genes.png', xlabel='Gene expression LFC, Tcon ÷ Treg')

    print("Bar p-values")
    for i, j in itertools.combinations(list(dict_for_pval.keys()), 2):
        (i_succes, i_tot), (j_success, j_tot) = dict_for_pval[i], dict_for_pval[j]
        pval = proportion_significance_test(i_succes, j_success, i_tot, j_tot)
        print(i, j, pval)

    ax.bar(names+["All genes"], sig_props+[(all_pvals<.05).mean()])
    ax.set_ylabel("Percent of genes with p < .05")
    fig.savefig('./plots/gene_cdf/pvals_bar.png')
    

import scipy
def binding_heatmap(ancdict, partners):
    sns.set(font_scale=1)
    for anc1_name, anc2_name in itertools.combinations(list(ancdict.keys()), 2):
        ancset1 = ancdict[anc1_name]
        ancset2 = ancdict[anc2_name]
        enrichment_df = pd.DataFrame()
        pval_df = pd.DataFrame()
        enrichments = []
        pvals = []
        for partner_name in list(partners.keys()):
            partner = partners[partner_name]


            partner_in_ancset1 = list(map(int, list(zip(*ancset1.intersect(partner, c=True)))[3]))
            partner_in_ancset2 = list(map(int, list(zip(*ancset2.intersect(partner, c=True)))[3]))

            p_in_ancset1 = np.mean(partner_in_ancset1)
            p_in_ancset2 = np.mean(partner_in_ancset2)
            if p_in_ancset2 != 0:
                enrichment = p_in_ancset1 / p_in_ancset2
                enrichments.append(enrichment)
                pval = scipy.stats.mannwhitneyu(partner_in_ancset1, partner_in_ancset2)
                pvals.append(pval[1])
            else:
                raise Exception
        enrichment_df[f'{anc1_name} / {anc2_name}'] = enrichments
        pval_df[f'{anc1_name} / {anc2_name}'] = pvals
        enrichment_df.index = list(partners.keys())
        pval_df.index = list(partners.keys())

        
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.heatmap(np.log2(enrichment_df),
           vmin=-1, vmax = 1, cmap=cm.bwr, annot=False, ax=ax)
        ax.figure.subplots_adjust(left = 0.3) 
        fig.savefig(f'./plots/binding_heatmap/{anc1_name}_vs_{anc2_name}_LFC.png')

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.heatmap(pval_df,
           vmin=0, vmax = .1, cmap=cm.bwr, annot=False, ax=ax)
        ax.figure.subplots_adjust(left = 0.3) 
        fig.savefig(f'./plots/binding_heatmap/{anc1_name}_vs_{anc2_name}_pval.png')


def one_sample_proportion_significance_test(p1_success, p1_tot, p2_success):
    zscore = (p1_success-p2_success)/(np.sqrt(p2_success*(1-p2_success))/p1_tot)
    pval = 2*(1-scipy.stats.norm.cdf(np.abs(zscore), 0, 1))    
    return pval

def proportion_significance_test(p1_success, p2_success, p1_tot, p2_tot):
    phat = (p1_success+p2_success)/(p1_tot + p2_tot)
    p1_prob = p1_success / p1_tot 
    p2_prob = p2_success / p2_tot
    zscore = (p1_prob-p2_prob)/np.sqrt(phat*(1-phat)*(1/p1_tot + 1/p2_tot))
    # pval = 2*(1-scipy.stats.norm.cdf(np.abs(zscore), 0, 1))    
    pval = 2*(1-scipy.stats.norm.cdf(np.abs(zscore), 0, 1))    
    return pval


def hicdc_vs_must(hicdc_only, Treg_hicdc, Tcon_hicdc, frame):
    for i in range(20):
        ind = np.random.choice(len(hicdc_only))
        place = hicdc_only[ind]

        chrom, start, end = place[:3]

        start, end = map(int, [start, end])
        anc = pbt.BedTool(f'{chrom} {start} {end}', from_string=True)
        p = Treg_hicdc.intersect(anc, u=True)
        q = Tcon_hicdc.intersect(anc, u=True)

        if len(p) > 0:
            chrom, start = p[0][:2]
            end = p[0][4]
        if len(q) > 0:
            chrom, start = q[0][:2]
            end = q[0][4]
        start, end = map(int, [start, end])
        start = start-20*5000
        end = end+20*5000
    
        frame *= Feature(depth_ratio=1)
        fig = frame.plot(chrom, start, end)

        fig.savefig(f'./plots/hicdc_vs_mustache/{i}.png')    


