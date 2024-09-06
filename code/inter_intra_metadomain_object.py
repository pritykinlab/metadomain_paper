from aux_functions import *
from make_figure4 import *

def trans(x):
    return list(zip(*x))

class MetadomainDataset:
    def __init__(self, all_intra_megaloops, all_inter_megaloops, all_ind_to_region, 
                        parsed_chroms, superenhancer_dict, bbi_dict, ind_to_gene, my_treg_comp, 
                        inter_and_intra_connections_treg, inter_and_intra_connections_tcon, cond='all',
                        n_initial_clusters = 15, random_state=2, intra_co=20, inter_co=20):
        self.intra_co = intra_co
        self.inter_co = inter_co
        self.inter_and_intra_connections_treg = inter_and_intra_connections_treg
        self.inter_and_intra_connections_tcon = inter_and_intra_connections_tcon
        self.my_treg_comp = my_treg_comp
        self.ind_to_gene = ind_to_gene
        self.bbi_dict = bbi_dict
        self.all_ind_to_region = all_ind_to_region
        self.parsed_chroms = parsed_chroms
        self.superenhancer_dict = superenhancer_dict
        self.all_intra_megaloops = all_intra_megaloops
        self.all_inter_megaloops = all_inter_megaloops
        self.all_megaloops = all_intra_megaloops + all_inter_megaloops
        self.cond = cond
        self.connect_dict = {
            'all' : self.all_megaloops > 0,
            'tcon' : self.inter_and_intra_connections_tcon > 0,
            'treg' : self.inter_and_intra_connections_treg > 0,
            # 'all_with_diag' : inter_and_intra_connections + np.diag(np.ones(len(inter_and_intra_connections))),
            # 'all_with_loops' : inter_and_intra_connections + all_short_range_loops,
        }

        self.binsoi = self.subset_for_clustering()
        self.submat_for_clustering = self.connect_dict[self.cond][self.binsoi, :][:, self.binsoi]
        # inter_and_intra_connections_tcon = all_intra_tcon_megaloops + all_inter_tcon_megaloops
        # inter_and_intra_connections_treg = all_intra_treg_megaloops + all_inter_treg_megaloops

        (out) = make_odict_and_clustdict_v2(self.connect_dict, self.parsed_chroms, self.all_ind_to_region,
                                            self.binsoi, n_clusters=n_initial_clusters, method='ward', n_pcs=20,
                                            random_state=random_state)
        clustdict, odict, le, goodinds = (out)

        self.clustdict = clustdict
        self.odict = odict
        self.le = le
        self.goodinds = arr(goodinds)
        self.chromlist = [self.all_ind_to_region[x] for x in self.goodinds]

        self.SE_count_dict = self.make_SE_dict(self.goodinds)
        self.bw_val_dict = self.make_bw_dict(self.goodinds)
        
        # Unordered submat
        self.focal_submat = self.get_focal_submat()
        (self.inds_with_enough_megaloops, self.cluster_to_subset_for_further_clustering, 
         self.n_megaloop_by_cluster_df, self.n_megaloops_within_cluster, 
         self.n_different_chroms_within_cluster, self.n_different_chroms_within_cluster_df) = self.reorder_and_subset_clusters()

        self.merged_clustdict, self.clusters_to_merge = self.merge_clusters()
        self.merged_inds_to_subset = [x for x  in self.cluster_to_subset_for_further_clustering if np.sum(self.merged_clustdict[self.cond]==x) > 0]
        self.merged_focal_submat = self.plot_inter_hub_focal_contacts(plot_ets1_se=False)
        
    def subset_for_clustering(self):
        intra_binsoi = set(np.where(self.all_intra_megaloops.sum(axis=1)>self.intra_co)[0])
        inter_binsoi = set(np.where(self.all_inter_megaloops.sum(axis=1)>self.inter_co)[0])
        binsoi = sorted(list(intra_binsoi.union(inter_binsoi)))
        return binsoi

    def make_SE_dict(self, goodinds):
        SE_count_dict = {}
        goodind_regions = pbt.BedTool([self.all_ind_to_region[i] for i in goodinds])
        for name, bedtool in self.superenhancer_dict.items():
            z = goodind_regions.intersect(bedtool, c=True)
            counts = get_col(z, -1).astype(float)
            SE_count_dict[name] = counts
        return SE_count_dict
    
    def make_bw_dict(self, goodinds):
        bw_val_dict = {}
        goodind_regions = pbt.BedTool([self.all_ind_to_region[i] for i in goodinds])
        for name, path in self.bbi_dict.items():
            bw = bbi.open(path)
            chrom, s, e = trans(add_chr_to_bedtool(goodind_regions))
            z = bw.stackup(chrom, s, e, bins=1)
            bw_val_dict[name] = z
        return bw_val_dict

    def plot_focal_contacts(self, ):
        fig, ax, focal_submat = plot_focal_contacts_final(self.connect_dict[self.cond],
                                                    self.all_ind_to_region, self.goodinds, self.odict[self.cond], 
                                                    self.clustdict[self.cond], self.le, self.my_treg_comp[self.goodinds], self.SE_count_dict, 
                                                    self.bw_val_dict, )
        return fig, ax, focal_submat

    def get_focal_submat(self, ):
        fig, ax, focal_submat = plot_focal_contacts_final(self.connect_dict[self.cond],
                                                    self.all_ind_to_region, self.goodinds, self.odict[self.cond], 
                                                    self.clustdict[self.cond], self.le,
                                                     self.my_treg_comp[self.goodinds], self.SE_count_dict, 
                                                    self.bw_val_dict, )
        fig.savefig('./plots/INTER_MEGALOOPS/focal_contacts.pdf', bbox_inches='tight')
        del fig, 
        return focal_submat


    def reorder_and_subset_clusters(self, ):
        ## Get chroms 
        chromlabels = arr(self.chromlist)[:, 0]
        o = self.odict[self.cond]
        ## Get chroms in order of clustering
        chrom_cluster_df = pd.DataFrame(chromlabels[o], 
                                                columns=['chrom'], 
                                                index=self.clustdict[self.cond][o]).reset_index().rename(
                                                columns={'index':'cluster'})
        ## Get number_of_megaloops_by_cluster
        chrom_cluster_df = chrom_cluster_df.value_counts().unstack()
        n_megaloop_by_cluster_df = pd.DataFrame(self.focal_submat).abs()
        n_megaloop_by_cluster_df['cluster'] = self.clustdict[self.cond]
        n_megaloop_by_cluster_df = n_megaloop_by_cluster_df.fillna(0).groupby('cluster').mean().T
        n_megaloop_by_cluster_df['cluster'] = self.clustdict[self.cond]
        n_megaloop_by_cluster_df = n_megaloop_by_cluster_df.groupby('cluster').mean()

        ## Get clusters to subset
        n_megaloops_within_cluster = pd.Series(np.diag(n_megaloop_by_cluster_df), n_megaloop_by_cluster_df.index)
        n_different_chroms_within_cluster_df = chrom_cluster_df.div(chrom_cluster_df.sum(axis=1), axis=0)
        n_different_chroms_within_cluster = n_different_chroms_within_cluster_df.max(axis=1)
        clusters_with_many_chroms = chrom_cluster_df[n_different_chroms_within_cluster < .3]
        inds_with_many_chroms = clusters_with_many_chroms.index
        inds_with_enough_megaloops = n_megaloops_within_cluster.index[n_megaloop_by_cluster_df.mean(axis=1) > 0.018]
        inds_to_subset_for_clustering = inds_with_many_chroms.intersection(inds_with_enough_megaloops)
        return (inds_with_enough_megaloops, inds_to_subset_for_clustering, n_megaloop_by_cluster_df,
                n_megaloops_within_cluster, n_different_chroms_within_cluster, n_different_chroms_within_cluster_df)

    def merge_clusters(self):
        merged_clustdict = deepcopy(self.clustdict)
        cluster_cormat = self.n_megaloop_by_cluster_df.corr()

        cluster_cormat = cluster_cormat.loc[self.cluster_to_subset_for_further_clustering, self.cluster_to_subset_for_further_clustering]   
        _, clusters_to_merge, _ = make_order_and_cluster_custom(cluster_cormat > .5, n_clusters=3)

        for u in np.unique(clusters_to_merge):
            indsoi = clusters_to_merge == u
            clusterlist = cluster_cormat.index[indsoi]
            if len(clusterlist) > 1:
                base_cluster = clusterlist[0]
                for i in clusterlist[1:]:
                    merged_clustdict[self.cond][self.clustdict[self.cond] == i] = base_cluster

        return merged_clustdict, clusters_to_merge

    def plot_inter_hub_focal_contacts(self, plot_ets1_se=True):
        genes_to_plot = ['Pdlim7', 'Lta', 'Jund', 'Cd4', 'Irf3', 'Cd37', 'Notch1', 'Lck', 
                        'Klf2', 'Lag3', 'Cd27', 'Coro1b', 'Foxa3', 'Ctcf', 'Sema4a', 'Entpd5', 
                        'Itgb4', 'Il17c', 'Cxcr5', 'Id3', 'Il4ra', 'Itga10', 'Fyn',
                        'Ctla4', 'Chd7', 'Bcl2',  'Tgfbr2',  'Tgfbr3', 'Il10', 'Lrrc56', 'Elf1', 'Cxcr4',  
                        'Btg1',  'Irf6', 'Fasl', 'Dusp10', 'Ptprc', 'Themis', 'Tox', 'Irf2',  'Wwox', 'Il7r',
                        'Runx1', 'Arl4c', 'Itk', 'Inpp4b', 'Foxp1', 
                        'Ets1', 'Lef1', 'Ikzf1', 'Ikzf4', 'Bach2',   
                        'Ikzf2', 'Il22', 'Izumo1r', 'Hivep2', 'Satb1',
                        'Stt3b', 'Top2a', 'Sgk1', 'Igf1r', 'Pde3b'
        ]
        cond = 'all'
        subdict = subset_dict(self.bw_val_dict, [ 'Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1', 'Treg Stat5',
                                            'Treg H3K27me3', ][::-1])

        fig, ax, merged_focal_submat, _ = plot_focal_contacts_submat_final(
                                                    self.inter_and_intra_connections_treg,
                                                    self.inter_and_intra_connections_tcon,
                                                    self.all_ind_to_region, self.goodinds, 
                                                    # merged_o, 
                                                    self.odict[self.cond], 
                                                    self.clustdict[self.cond],
                                                    self.merged_clustdict[self.cond], self.le, 
                                                    self.my_treg_comp[self.goodinds], 
                                                    self.SE_count_dict, 
                                                    subdict, self.merged_inds_to_subset, 
                                                    genes_to_plot = genes_to_plot, 
                                                    ind_to_gene = self.ind_to_gene,
                                                    plot_ets1_se=plot_ets1_se
                                                    )	
        fig.savefig('./plots/FINAL_MEGALOOPS/inter_hub.svg', bbox_inches='tight')
        return merged_focal_submat	


    def plot_inter_hub_focal_contacts_baseline(self):
        genes_to_plot = ['Pdlim7', 'Lta', 'Jund', 'Cd4', 'Irf3', 'Cd37', 'Notch1', 'Lck', 
                        'Klf2', 'Lag3', 'Cd27', 'Coro1b', 'Foxa3', 'Ctcf', 'Sema4a', 'Entpd5', 
                        'Itgb4', 'Il17c', 'Cxcr5', 'Id3', 'Il4ra', 'Itga10', 'Fyn',
                        'Ctla4', 'Chd7', 'Bcl2',  'Tgfbr2',  'Tgfbr3', 'Il10', 'Lrrc56', 'Elf1', 'Cxcr4',  
                        'Btg1',  'Irf6', 'Fasl', 'Dusp10', 'Ptprc', 'Themis', 'Tox', 'Irf2',  'Wwox', 'Il7r',
                        'Runx1', 'Arl4c', 'Itk', 'Inpp4b', 'Foxp1', 
                        'Ets1', 'Lef1', 'Ikzf1', 'Ikzf4', 'Bach2',   
                        'Ikzf2', 'Il22', 'Izumo1r', 'Hivep2', 'Satb1',
                        'Stt3b', 'Top2a', 'Sgk1', 'Igf1r', 'Pde3b'
        ]
        cond = 'all'
        subdict = subset_dict(self.bw_val_dict, [ 'Treg H3K27ac', 'Treg H3K4me3', 'Treg H3K4me1', 'Treg Stat5',
                                            'Treg H3K27me3', ][::-1])

        fig, ax, merged_focal_submat, _ = plot_focal_contacts_submat_final(
                                                    self.connect_dict[self.cond],
                                                    np.zeros_like(self.connect_dict[self.cond]),
                                                    self.all_ind_to_region, self.goodinds, 
                                                    # merged_o, 
                                                    self.odict[self.cond], 
                                                    self.merged_clustdict[self.cond], self.le, 
                                                    self.my_treg_comp[self.goodinds], 
                                                    self.SE_count_dict, 
                                                    subdict, self.merged_inds_to_subset, 
                                                    genes_to_plot = genes_to_plot, 
                                                    ind_to_gene = self.ind_to_gene
                                                    )	
        fig.savefig('./plots/FINAL_MEGALOOPS/inter_hub.svg', bbox_inches='tight')
        return merged_focal_submat	

    def plot_reorder_inter_hub_focal_contacts(self):
        genes_to_plot = ['Pdlim7', 'Lta', 'Jund', 'Cd4', 'Irf3', 'Cd37', 'Notch1', 'Lck', 
                        'Klf2', 'Lag3', 'Cd27', 'Coro1b', 'Foxa3', 'Ctcf', 'Sema4a', 'Entpd5', 
                        'Itgb4', 'Il17c', 'Cxcr5', 'Id3', 'Il4ra', 'Itga10', 'Fyn',
                        'Ctla4', 'Chd7', 'Bcl2',  'Tgfbr2',  'Tgfbr3', 'Il10', 'Lrrc56', 'Elf1', 'Cxcr4',  
                        'Btg1',  'Irf6', 'Fasl', 'Dusp10', 'Ptprc', 'Themis', 'Tox', 'Irf2',  'Wwox', 'Il7r',
                        'Runx1', 'Arl4c', 'Itk', 'Inpp4b', 'Foxp1', 'Ets1', 'Lef1', 'Ikzf1', 'Ikzf4', 'Bach2',   
                        'Ikzf2', 'Il22', 'Izumo1r', 'Hivep2', 'Satb1',
                        'Stt3b', 'Top2a', 'Sgk1', 'Igf1r', 'Pde3b',
                        'Rxra', 'Rbfox3', 'Lrrc16b', 'Sox8'
        ]
        cond = 'all'
        subdict = subset_dict(self.bw_val_dict, ['Treg Stat5', 'Treg H3K4me3', 'Treg H3K4me1',  'Treg H3K27ac', 
                                            'Treg H3K27me3', 'Treg CTCF', 'Treg Smc1a', ])

        subinds = np.isin(selc.merged_clustdict[self.cond], merged_inds_to_subset)
        merged_o, _1, _ = make_order_and_cluster_custom(self.focal_submat[self.goodinds, :][:, self.goodinds], n_clusters=3, 
                                                        metric='euclidean', method='ward')
        fig, ax, merged_focal_submat, _ = plot_focal_contacts_submat_final(
                                                    self.inter_and_intra_connections_treg,
                                                    self.inter_and_intra_connections_tcon,
                                                    self.all_ind_to_region, self.goodinds, 
                                                    merged_o, 
                                                    # self.odict[self.cond], 
                                                    self.merged_clustdict[self.cond], self.le, 
                                                    self.my_treg_comp[self.goodinds], 
                                                    self.SE_count_dict, 
                                                    subdict, self.merged_inds_to_subset, 
                                                    genes_to_plot = genes_to_plot, 
                                                    ind_to_gene = self.ind_to_gene
                                                    )	
        return merged_focal_submat	


def aggregate_matrix_by_clusters_and_inds_to_subset(matrix_to_cluster, inds_to_subset, clusts, o, ignore_nan=True):
    n = len(np.unique(inds_to_subset))
    cluster_oe = np.zeros((n, n))
    for cl1 in range(n):
        for cl2 in range(n):
            indsoi1 = clusts[o] == inds_to_subset[cl1]
            indsoi2 = clusts[o] == inds_to_subset[cl2]
            if ignore_nan:
                avg_oe = np.nanmean(matrix_to_cluster[indsoi1, :][:, indsoi2])
            else:
                m = matrix_to_cluster[indsoi1, :][:, indsoi2].copy()
                m[np.isnan(m)] = 0
                avg_oe = np.nanmean(m)
            cluster_oe[cl1, cl2] = avg_oe

    cluster_oe = pd.DataFrame(cluster_oe)
    cluster_oe.index = inds_to_subset
    cluster_oe.columns = inds_to_subset    
    return cluster_oe


def aggregate_matrix_by_clusters(matrix_to_cluster, clusts, ignore_nan=True):
    n = len(np.unique(clusts))
    cluster_oe = np.zeros((n, n))
    for cl1 in range(n):
        for cl2 in range(n):
            indsoi1 = clusts == cl1
            indsoi2 = clusts == cl2
            if ignore_nan:
                avg_oe = np.nanmean(matrix_to_cluster[indsoi1, :][:, indsoi2])
            else:
                m = matrix_to_cluster[indsoi1, :][:, indsoi2].copy()
                m[np.isnan(m)] = 0
                avg_oe = np.nanmean(m)
            cluster_oe[cl1, cl2] = avg_oe

    cluster_oe = pd.DataFrame(cluster_oe)
    cluster_oe.index = range(n)
    cluster_oe.columns = range(n)    
    return cluster_oe


def get_dfs_with_cluster_megaloopcounts(inter_and_intra_connections, sep_oe_mat_tcon, sep_oe_mat_treg, megaloop_dataset):
    o = megaloop_dataset.odict[megaloop_dataset.cond]
    goodinds = megaloop_dataset.goodinds
    oe_tcon_here = sep_oe_mat_tcon[goodinds, :][o][:, goodinds][:, o]
    oe_treg_here = sep_oe_mat_treg[goodinds, :][o][:, goodinds][:, o]
    inds_to_subset = megaloop_dataset.merged_inds_to_subset
    cond = megaloop_dataset.cond
    clustdict = megaloop_dataset.merged_clustdict
    # megaloops_tcon_in_cluster = megaloop_dataset.inter_and_intra_connections_tcon[goodinds, :][o][:, goodinds][:, o]
    # megaloops_treg_in_cluster = megaloop_dataset.inter_and_intra_connections_treg[goodinds, :][o][:, goodinds][:, o]

    # tmp = inter_and_intra_connections[goodinds, :][:, goodinds][o, :][:, o]

    # cluster_oe_tcon = aggregate_matrix_by_clusters_and_inds_to_subset(oe_tcon_here, inds_to_subset, clustdict[cond], o)
    # cluster_oe_treg = aggregate_matrix_by_clusters_and_inds_to_subset(oe_treg_here, inds_to_subset, clustdict[cond], o)

    # cluster_megaloop_tcon = aggregate_matrix_by_clusters_and_inds_to_subset(megaloops_tcon_in_cluster, inds_to_subset, 
    #                                                     clustdict[cond], o, ignore_nan = False)
    # cluster_megaloop_treg = aggregate_matrix_by_clusters_and_inds_to_subset(megaloops_treg_in_cluster, inds_to_subset, 
    #                                                     clustdict[cond], o, ignore_nan = False)	
    tcon_megaloop_df = pd.DataFrame(megaloop_dataset.inter_and_intra_connections_tcon[:, goodinds][:, o], columns=goodinds[o])
    treg_megaloop_df = pd.DataFrame(megaloop_dataset.inter_and_intra_connections_treg[:, goodinds][:, o], columns=goodinds[o])
    
    megaloop_cluster_count_treg = pd.DataFrame()
    megaloop_cluster_count_tcon = pd.DataFrame()
    for cluster in inds_to_subset:
        indsoi = clustdict[cond][o] == cluster
        if indsoi.sum() == 0:
            continue
        treg_megaloops_in_cluster = treg_megaloop_df.loc[:, indsoi]
        tcon_megaloops_in_cluster = tcon_megaloop_df.loc[:, indsoi]
        
        x = treg_megaloops_in_cluster.sum(axis=1)
        y = tcon_megaloops_in_cluster.sum(axis=1)
        megaloop_cluster_count_treg[cluster] = x
        megaloop_cluster_count_tcon[cluster] = y

    return megaloop_cluster_count_treg, megaloop_cluster_count_tcon


def get_dfs_with_cluster_hic_oe(inter_and_intra_connections, sep_oe_mat_tcon, sep_oe_mat_treg, megaloop_dataset):
    o = megaloop_dataset.odict[megaloop_dataset.cond]
    goodinds = megaloop_dataset.goodinds
    oe_tcon_here = sep_oe_mat_tcon[goodinds, :][o][:, goodinds][:, o]
    oe_treg_here = sep_oe_mat_treg[goodinds, :][o][:, goodinds][:, o]
    inds_to_subset = megaloop_dataset.merged_inds_to_subset
    cond = megaloop_dataset.cond
    clustdict = megaloop_dataset.merged_clustdict
    megaloops_tcon_in_cluster = megaloop_dataset.inter_and_intra_connections_tcon[goodinds, :][o][:, goodinds][:, o]
    megaloops_treg_in_cluster = megaloop_dataset.inter_and_intra_connections_treg[goodinds, :][o][:, goodinds][:, o]

    tmp = inter_and_intra_connections[goodinds, :][:, goodinds][o, :][:, o]


    cluster_oe_tcon = aggregate_matrix_by_clusters_and_inds_to_subset(oe_tcon_here, inds_to_subset, clustdict[cond], o)
    cluster_oe_treg = aggregate_matrix_by_clusters_and_inds_to_subset(oe_treg_here, inds_to_subset, clustdict[cond], o)

    cluster_megaloop_tcon = aggregate_matrix_by_clusters_and_inds_to_subset(megaloops_tcon_in_cluster, inds_to_subset, 
                                                        clustdict[cond], o, ignore_nan = False)
    cluster_megaloop_treg = aggregate_matrix_by_clusters_and_inds_to_subset(megaloops_treg_in_cluster, inds_to_subset, 
                                                        clustdict[cond], o, ignore_nan = False)	
    tcon_oe_df = pd.DataFrame(sep_oe_mat_tcon[:, goodinds][:, o], columns=goodinds[o])
    treg_oe_df = pd.DataFrame(sep_oe_mat_treg[:, goodinds][:, o], columns=goodinds[o])
    
    megaloop_cluster_oe_treg = pd.DataFrame()
    megaloop_cluster_oe_tcon = pd.DataFrame()
    for cluster in inds_to_subset:
        indsoi = clustdict[cond][o] == cluster
        if indsoi.sum() == 0:
            continue
        treg_oe_in_cluster = treg_oe_df.loc[:, indsoi]
        tcon_oe_in_cluster = tcon_oe_df.loc[:, indsoi]
        
        x = treg_oe_in_cluster.mean(axis=1)
        y = tcon_oe_in_cluster.mean(axis=1)
        megaloop_cluster_oe_treg[cluster] = x
        megaloop_cluster_oe_tcon[cluster] = y

    return megaloop_cluster_oe_treg, megaloop_cluster_oe_tcon

def get_dfs_with_differential_megaloops_to_clusters(inter_and_intra_connections, sep_oe_mat_tcon, sep_oe_mat_treg, megaloop_dataset):
    o = megaloop_dataset.odict[megaloop_dataset.cond]
    goodinds = megaloop_dataset.goodinds
    oe_tcon_here = sep_oe_mat_tcon[goodinds, :][o][:, goodinds][:, o]
    oe_treg_here = sep_oe_mat_treg[goodinds, :][o][:, goodinds][:, o]
    inds_to_subset = megaloop_dataset.merged_inds_to_subset
    cond = megaloop_dataset.cond
    clustdict = megaloop_dataset.merged_clustdict
    megaloops_tcon_in_cluster = megaloop_dataset.inter_and_intra_connections_tcon[goodinds, :][o][:, goodinds][:, o]
    megaloops_treg_in_cluster = megaloop_dataset.inter_and_intra_connections_treg[goodinds, :][o][:, goodinds][:, o]

    tmp = inter_and_intra_connections[goodinds, :][:, goodinds][o, :][:, o]


    cluster_oe_tcon = aggregate_matrix_by_clusters_and_inds_to_subset(oe_tcon_here, inds_to_subset, clustdict[cond], o)
    cluster_oe_treg = aggregate_matrix_by_clusters_and_inds_to_subset(oe_treg_here, inds_to_subset, clustdict[cond], o)

    cluster_megaloop_tcon = aggregate_matrix_by_clusters_and_inds_to_subset(megaloops_tcon_in_cluster, inds_to_subset, 
                                                        clustdict[cond], o, ignore_nan = False)
    cluster_megaloop_treg = aggregate_matrix_by_clusters_and_inds_to_subset(megaloops_treg_in_cluster, inds_to_subset, 
                                                        clustdict[cond], o, ignore_nan = False)	


    tcon_megaloop_df = pd.DataFrame(megaloop_dataset.inter_and_intra_connections_tcon[:, goodinds][:, o], columns=goodinds[o])
    treg_megaloop_df = pd.DataFrame(megaloop_dataset.inter_and_intra_connections_treg[:, goodinds][:, o], columns=goodinds[o])

    (megaloop_frequency_enrichment_df, 
        megaloop_frequency_pval_df) = get_cluster_specific_megaloop_enrichments_by_counting(
                                            treg_megaloop_df, tcon_megaloop_df, 
                                            inds_to_subset, clustdict, o)

    (megaloop_frequency_enrichment_df_permute, 
        megaloop_frequency_pval_df_permute) = get_cluster_specific_megaloop_enrichments_by_permuting(
                                            treg_megaloop_df, tcon_megaloop_df, 
                                            inds_to_subset, clustdict, o)
    return megaloop_frequency_enrichment_df, megaloop_frequency_pval_df, megaloop_frequency_enrichment_df_permute, megaloop_frequency_pval_df_permute

def get_dfs_with_differential_interactions_to_clusters(inter_and_intra_connections, inter_oe_mat_treg, inter_oe_mat_tcon, megaloop_dataset):
    o = megaloop_dataset.odict[megaloop_dataset.cond]
    goodinds = megaloop_dataset.goodinds
    
    treg_interaction_df = pd.DataFrame(inter_oe_mat_treg[:, goodinds][:, o], columns=goodinds[o])
    tcon_interaction_df = pd.DataFrame(inter_oe_mat_tcon[:, goodinds][:, o], columns=goodinds[o])
    # return treg_interaction_df, tcon_interaction_df
    inds_to_subset = megaloop_dataset.merged_inds_to_subset
    clustdict = megaloop_dataset.merged_clustdict

    (cluster_hic_enrichment_df, cluster_hic_pval_df) = get_cluster_specific_megaloop_enrichments_by_permuting(treg_interaction_df, 
                                                                                                            tcon_interaction_df, 
                                                                                                            inds_to_subset, 
                                                                                                            clustdict, 
                                                                                                            o)	
    return cluster_hic_enrichment_df, cluster_hic_pval_df
## Remove artifactual chromosome
# cluster_to_subset_for_further_clustering = cluster_to_subset_for_further_clustering.drop(22)		