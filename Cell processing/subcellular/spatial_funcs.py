from numpy.core.shape_base import vstack
import scanpy as sc
import scipy.io as sio
from starmap.sequencing import *

from sklearn.neighbors import NearestNeighbors
from scipy.spatial import Delaunay
from sklearn.cluster import AgglomerativeClustering as ac
import faiss

import numpy as np
import pandas as pd
import math
from scipy.stats import norm, binom
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
from scipy.ndimage import distance_transform_edt as dt
from skimage import measure

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

from copy import copy
import random
from tqdm import tqdm

# from reads assignment
def load_reads_2D(fpath, reads_file):
    S = loadmat(os.path.join(fpath, reads_file))
    bases = [str(i[0][0]) for i in S["merged_reads"]]
    points = S["merged_points"][:, :2]
    temp = np.zeros(points.shape)
    temp[:, 0] = np.round(points[:, 1]-1)
    temp[:, 1] = np.round(points[:, 0]-1)

    print(f"Number of reads: {len(bases)}")

    return bases, temp

def load_reads_3D(fpath, reads_file):
    S = loadmat(os.path.join(fpath, reads_file))
    bases = [str(i[0][0]) for i in S["merged_reads"]]
    points = S["merged_points"]
    temp = np.zeros(points.shape)
    temp[:, 0] = np.round(points[:, 2]-1)
    temp[:, 1] = np.round(points[:, 1]-1)
    temp[:, 2] = np.round(points[:, 0]-1)

    print(f"Number of reads: {len(bases)}")

    return bases, temp

def load_genes(fpath):
    genes2seq = {}
    seq2genes = {}
    with open(os.path.join(fpath, "genes.csv"), encoding='utf-8-sig') as f:
        for l in f:
            fields = l.rstrip().split(",")
            genes2seq[fields[0]] = "".join([str(s+1) for s in encode_SOLID(fields[1][::-1])])
            seq2genes[genes2seq[fields[0]]] = fields[0]
    return genes2seq, seq2genes

def points_by_gene(points, bases, gene):
    return points[np.argwhere(bases == gene)[:,0]]

def subset_points(pts, upper_bound):
    ### for points of 2d array; plotting purposes
    return pts[np.argwhere(np.logical_and(pts[:,0] < upper_bound, pts[:,1] < upper_bound))[:,0]]

def find_knn(gene_list, cell_list, points, bases, reads_to_cell, n_neigh=10 ,n_perm=0): # n_jobs, radius
    nn = NearestNeighbors(n_neighbors = n_neigh)
    dist = {gene: None for gene in gene_list}
    for gene in tqdm(gene_list):
        counts = {gene: 0 for gene in gene_list} # each [] contains a num
        permutation = {gene: [0] * n_perm for gene in gene_list} # each [] contains n_perm counts
        for cell in cell_list: # np.unique(reads_to_cell)[1:]: # 412 num per loop
            pts = np.argwhere(reads_to_cell == cell)[:,0]
            fitting_reads = np.argwhere(bases[pts] != gene)[:,0]
            query_reads = np.argwhere(bases[pts] == gene)[:,0]
            if query_reads.shape[0] < 5 or fitting_reads.shape[0] < 200: ### threshold # NOTE maybe remove if adding radius
                continue
            nn.fit(points[pts][fitting_reads])
            I = nn.kneighbors(points[pts][query_reads], return_distance = False) # NOTE if add radius -- add 'if I:'
            input_id = bases[pts][fitting_reads] # input id (into nn)
            obs = input_id[I]
            # count
            gene_id, ct = np.unique(obs, return_counts = True)
            for i, val in enumerate(gene_id):
                counts[val] += ct[i]

            # permutation
            if n_perm > 0:
                for p in range(n_perm):
                    np.random.shuffle(input_id) # shuffling happens in place
                    perm = input_id[I] ### TBD -- current shuffling: shuffle labels in the cell which != randomly assign label
                    ge_id, cts = np.unique(perm, return_counts = True)
                    for j, value in enumerate(ge_id):
                        permutation[value][p] += cts[j]

            dist[gene] = counts, permutation
        #     break
        # break
    return dist

def pval_dist(dist, gene_list, n_perm):
    '''Do not check whether NaN or None exists'''
    pVal = []
    for g in gene_list:
        pV = []
        for gg in gene_list:
            pV.append(np.count_nonzero(np.array(dist[g][1][gg]) >= dist[g][0][gg]) / n_perm)
        pVal.append(pV)
    pVal = np.array(pVal)
    np.fill_diagonal(pVal, 1)
    return pVal

def viz_nn(dist_matrix, n_perm, gene_list, sanityCheck=False, sCgenes=None, heatMap=False, clusterMap=False, row_colors=None, subCluster=False, pCut=0.05, n_sigNN=10):
    '''Could be used to just return pval matrix (when dist_matrix has None entries)
        row_colors should be in the same order as gene_list'''

    # sanity check
    if sanityCheck:
        p = sns.distplot(dist_matrix[sCgenes[0]][1][sCgenes[1]], color = 'g', kde = False)
        plt.axvline(dist_matrix[sCgenes[0]][0][sCgenes[1]])
        p.set(xlabel = '# times being kNN', ylabel = 'counts')

    # p-value matrix
    pVal = []
    print('No cells passed threshold for gene: ')
    num_fail = 0
    recorded_genes = []
    for g in gene_list:
        if dist_matrix[g] is None:
            print(g, ', ', end = '')
            num_fail += 1
            continue
        recorded_genes.append(g)
        pV = []
        for gg in gene_list:
            pV.append(np.count_nonzero(np.array(dist_matrix[g][1][gg]) >= dist_matrix[g][0][gg]) / n_perm)
        pVal.append(pV)
    print(num_fail, 'genes in total.')
    pVal = np.array(pVal)
    np.fill_diagonal(pVal, 1)

    # heatmap
    if heatMap:
        plt.figure(figsize = (15,15))
        sns.heatmap(pVal, cmap = 'YlGnBu', cbar_kws={'label': 'p values'})

    # clustermap
    if clusterMap:
        row_colors = np.array(row_colors)
        row_colors = row_colors[np.isin(gene_list, recorded_genes)]
        sns.clustermap(pVal, method='ward', cmap='mako', row_colors=row_colors, col_cluster=True, cbar_kws={'label': 'p values'})

    # subcluster
    if subCluster:
        sig_genes = []
        sig_pVal = []
        for rg in range(len(recorded_genes)):
            if np.count_nonzero(pVal[rg] < pCut) >= n_sigNN:
                sig_genes.append(recorded_genes[rg])
                sig_pVal.append(pVal[rg])
        plt.figure(figsize = (20,20))
        sns.clustermap(sig_pVal, method='ward', metric='euclidean', cmap='YlGnBu', col_cluster=True,
                        cbar_kws={'label': 'p values'}, yticklabels=sig_genes) #, row_colors = ['orange', 'yellow', 'red']

    return pVal, recorded_genes

def close_gene_pairs(pval_matrix, gene_list, n_pairs):
    thres = np.sort(pval_matrix, axis=None)[n_pairs-1]
    pair_loc = np.argwhere(pval_matrix <= thres)
    pairs = []
    for pair in pair_loc:
        if (pair[1], pair[0]) in pairs:
            continue
        else: pairs.append(tuple(pair))
    pairs = np.array(pairs)
    df_pairs = pd.DataFrame(zip(np.array(gene_list)[pairs[:,0]], np.array(gene_list)[pairs[:,1]], pval_matrix[pairs[:,0], pairs[:,1]]), columns=['Gene 1', 'Gene 2', 'P value'])
    df = df_pairs.sort_values(by=['P value'], ignore_index=True)
    return df

def close_neighbors(pval_matrix, genes_col, n_neigh=10): # TODO row_genes, col_genes
    close_neigh = []
    for p in pval_matrix:
        p_gene = sorted(list(zip(p, genes_col)))
        close = p_gene[:n_neigh]
        close_ = [c[1] for c in close]
        close_neigh.append(close_)
    return np.array(close_neigh)

def intersect_gene_pairs(df1, df2):
    '''pairs1/2 are dataframes with columns=['Gene1', 'Gene2', 'P value']
        gene pairs are sorted alphabetically to avoid (A, B) and (B, A) recognized as different'''
    intersect = []
    pair1 = list(zip(df1['Gene 1'], df1['Gene 2']))
    pair2 = list(zip(df2['Gene 1'], df2['Gene 2']))
    for p in pair1:
        if p in pair2:
            intersect.append(p)
        elif [p[1], p[0]] in pair2:
            intersect.append(p)
    return len(intersect), intersect

def close_gene_table(pval_matrix, genes_row, genes_col, n_neigh=10):
    close_neigh = close_neighbors(pval_matrix, genes_col, n_neigh)
    # col_name = ['Gene']
    col_name = ['neighbor-'+str(i) for i in range(n_neigh)]
    df = pd.DataFrame(close_neigh, columns=col_name)
    df.index = genes_row
    return df

def jaccard_graph(genes, metric, threshold, node_color='#210070'):
    # table of n_gene x n_gene of jaccard index values
    jaccard = np.zeros((len(genes),len(genes)))
    for i in range(len(genes)-1):
        for j in range(i+1,len(genes)):
            pair_jaccard = len(np.intersect1d(metric[i,:], metric[j,:]))/len(np.union1d(metric[i,:], metric[j,:]))
            jaccard[i,j] = pair_jaccard
            jaccard[j,i] = pair_jaccard

    # create a graph
    G = nx.Graph()

    # add node
    G.add_nodes_from(genes)

    # add edge
    for i in range(len(genes)-1):
        for j in range(i+1,len(genes)):
            if jaccard[i,j] > threshold: ## drawing all edges too messy -> filter out low similarity ones (weight=0 -- still a lot; weight<0.25 - two shared neighbors)
                G.add_edge(genes[i], genes[j], weight=jaccard[i,j])

    print('n_nodes: ', G.number_of_nodes(), '\nn_edges: ', G.number_of_edges())

    # set edge_width to weight (Jaccard)
    edge_width = [G.get_edge_data(u,v)['weight'] for u,v in G.edges]
    edge_width = np.array(edge_width)
    adj_edge_width = (edge_width * 10)**2

    # set node size to: sum of edge weights of connected edges
    node_size = []
    for gene in genes:
        node_size.append(sum([G.get_edge_data(u,v)['weight'] for u,v in G.edges(gene)]))
    node_size = np.array(node_size)
    adj_node_size = (node_size * 10)**2.4

    plt.figure(figsize=(12,12))
    pos = nx.kamada_kawai_layout(G) # spring-force-directed layout
    nx.draw_networkx_nodes(G, pos, node_size=adj_node_size, node_color=node_color, alpha=0.7)
    nx.draw_networkx_edges(G, pos, width=adj_edge_width, edge_color="#D9D9D9", alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8, bbox={"ec": "k", "fc": "white", "alpha": 0.7})
    return

def get_pts_bs_rtc(points, bases, genes, cells):
    pts = points[np.isin(bases, genes)]
    bs = bases[np.isin(bases, genes)]
    rtc = cells[pts[:,0], pts[:,1]]
    return pts, bs, rtc

def get_pts_bs_rtc_3D(points, bases, genes, cells):
    pts = points[np.isin(bases, genes)]
    bs = bases[np.isin(bases, genes)]
    rtc = cells[pts[:,0], pts[:,1], pts[:,2]]
    return pts, bs, rtc

def dl_nn(gene_list, cell_list, points, bases, reads_to_cell, n_perm=0):
    # nn = NearestNeighbors(n_neighbors = n_neigh)
    dist = {gene: None for gene in gene_list}
    for gene in gene_list:
        counts = {gene: 0 for gene in gene_list} # each [] contains a num
        permutation = {gene: [0] * n_perm for gene in gene_list} # each [] contains n_perm counts
        for cell in cell_list: # np.unique(reads_to_cell)[1:]: # 412 num per loop
            pts = np.argwhere(reads_to_cell == cell)[:,0]
            # fitting_reads = np.argwhere(bases[pts] != gene)[:,0]
            query_reads = np.argwhere(bases[pts] == gene)[:,0]
            if query_reads.shape[0] < 5 or fitting_reads.shape[0] < 200: ### threshold
                continue
            tri = Delaunay(points[pts]) # instead of only fitting reads because need to construct triangles with query reads ### post-processed in pval_DLdist
            # I = nn.kneighbors(points[pts][query_reads], return_distance = False)
            # input_id = bases[pts][fitting_reads] # input id (into nn)
            # obs = input_id[I]
            indptr, indices = tri.vertex_neighbor_vertices
            nn = [indices[indptr[k]:indptr[k+1]] for k in query_reads] # query_reads: index wrt points[pts]
            input_id = bases[pts]
            obs = list()
            for n in nn:
                obs.extend(list(input_id[n]))
            # count
            gene_id, ct = np.unique(obs, return_counts = True)
            for i, val in enumerate(gene_id):
                counts[val] += ct[i]

            # permutation
            if n_perm > 0:
                for p in range(n_perm):
                    np.random.shuffle(input_id) # shuffling happens in place
                    # perm = bases[pts][I] ### TBD -- current shuffling: shuffle labels in the cell which != randomly assign label
                    perm = list()
                    for n in nn:
                        perm.extend(list(input_id[n]))
                    ge_id, cts = np.unique(perm, return_counts = True)
                    for j, value in enumerate(ge_id):
                        permutation[value][p] += cts[j]

            dist[gene] = counts, permutation
    return dist

def draw_subset(points, bases, genes, cmap, seg_mask, upper_bound, dot_size=10):
    '''Draw part of the plot'''
    if len(genes) != len(cmap):
        raise ValueError('check ur input kiddo')
    # draw the mask
    plt.figure(figsize=(20,20))
    plt.imshow(seg_mask[:upper_bound, :upper_bound] > 0, cmap='gray')
    # get subset of reads and draw
    for i in range(len(genes)):
        curr_reads = subset_points(points_by_gene(points, bases, genes[i]), upper_bound)
        plt.scatter(curr_reads[:,0], curr_reads[:,1], s=dot_size, c=cmap[i])

def normalize_layer_by_cell_total(adata, layers):
    '''Normalize by cell total -- divide entry in layers by total_reads in cell (not in compartment)'''
    adata_new = copy(adata)
    for layer in layers:
        adata_new.layers[layer+'_norm'] = adata_new.layers[layer] / np.sum(adata_new.X, axis=1)[:,None]
    return adata_new

def adata_vstack_layers(adata, layers, attribute):
    '''vstack adata.X, keep adata.var unchanged, add attribute column to adata.obs'''
    X_new = np.vstack(tuple([adata.layers[layer] for layer in layers]))
    new_obs = []
    for layer in layers:
        obs = copy(adata.obs)
        obs[attribute] = layer
        new_obs.append(obs)
    obs_new = pd.concat(new_obs)
    adata_new = sc.AnnData(X=X_new, var=adata.var, obs=obs_new)
    adata_new.obs_names_make_unique()
    return adata_new

def adata_vstack_layers_norm(adata, layers, attribute):
    '''with normalize by cell total'''
    adata_new = normalize_layer_by_cell_total(adata, layers)
    layers_norm = [layer+'_norm' for layer in layers]
    adata_new = adata_vstack_layers(adata_new, layers_norm, attribute)
    return adata_new

def adata_hstack_layers(adata, layers):
    X_h = np.hstack(tuple(adata.layers[layer] for layer in layers))
    var_h = pd.concat([adata.var for i in range(len(layers))])
    adata_h = sc.AnnData(X=X_h, var=var_h, obs=adata.obs)
    adata_h.obs_names_make_unique()
    return adata_h

def separate_adata_by_sample(adata, sample_list, adata_layers=None):
    '''for separating AnnData that has all samples (with layers)
        where var is kept unchanged while obs & layers are separated'''
    X = adata.X[np.isin(adata.obs['sample'], sample_list)]
    obs = adata.obs.loc[np.isin(adata.obs['sample'], sample_list)]
    adata_new = sc.AnnData(X=X, obs=obs, var=adata.var)
    if adata_layers:
        for layer in adata_layers:
            adata_new.layers[layer] = adata.layers[layer][np.isin(adata.obs['sample'], sample_list)]
    return adata_new

def merge_adata_by_sample(adata, sample_list, adata_layers=None):
    '''adata being a dict of adata's to be merged, assuming same vars'''
    X = np.concatenate([adata[sample].X for sample in sample_list])
    obs = pd.concat([adata[sample].obs for sample in sample_list])
    adata_new = sc.AnnData(X=X, obs=obs, var=adata[sample_list[0]].var)
    if adata_layers:
        for layer in adata_layers:
            adata_new.layers[layer] = np.concatenate([adata[sample].layers[layer] for sample in sample_list])
    return adata_new

def adata_remove_genes(adata, genes_to_remove, layers=None):
    var = adata.var.loc[np.isin(adata.var.index, genes_to_remove, invert=True)]
    X = adata.X[:, np.isin(adata.var.index, genes_to_remove, invert=True)]
    adata_new = sc.AnnData(X=X, obs=adata.obs, var=var)
    if layers:
        for layer in layers:
            adata_new.layers[layer] = adata.layers[layer][:, np.isin(adata.var.index, genes_to_remove, invert=True)]
    return adata_new

def adata_remove_cells(adata, cell_label_category, cell_label, layers=None):
    '''remove cells with given label and adjust obs accordingly'''
    obs = adata.obs.loc[np.isin(adata.obs[cell_label_category], cell_label, invert=True)]
    X = adata.X[np.isin(adata.obs[cell_label_category], cell_label, invert=True), :]
    adata_new = sc.AnnData(X=X, obs=obs, var=adata.var)
    if layers:
        for layer in layers:
            adata_new.layers[layer] = adata.layers[layer][np.isin(adata.obs[cell_label_category], cell_label, invert=True), :]
    return adata_new

def adata_integration(adata_multiple, sample, int_col=None):
    '''simply vertically stacking two adata's Xs;
        original adata need to have sample in obs'''
    obs_int = []
    for _, adata in enumerate(adata_multiple):
        obs = copy(adata.obs)
        if int_col:
            obs['integration'] = obs[int_col].astype(str)+'_'+obs[sample].astype(str)
        obs_int.append(obs)
    adata_int = sc.AnnData(X=np.vstack(tuple(adata.X for adata in adata_multiple)),
                            var=pd.DataFrame(adata_multiple[0].var.index),
                            obs=pd.concat(obs_int))
    adata_int.obs_names_make_unique()
    return adata_int

def pp_pca_umap(adata, pp=True, layer=None, labels=None, pca_plot=False, pca_var=False, pca_loading=False, tsne=False, umap=False):
    # standard pre-processing
    if pp:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, zero_center=True)
    # pca
    sc.pp.pca(adata)
    if pca_plot:
        sc.pl.pca(adata, color=labels)
    if pca_var:
        sc.pl.pca_variance_ratio(adata)
    if pca_loading:
        sc.pl.pca_loadings(adata)
    # tsne
    if tsne:
        sc.tl.tsne(adata)
        sc.pl.tsne(adata, color=labels)
    # umap
    if umap:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=labels)

    return

def gene_by_compartment_from_vstack(adata, adata_v, layers):
    '''turn v-stacked adata.X into gene by compartment (layer) dataframe'''
    sum_cells = []
    for i in range(len(layers)):
        sum_cells.append(np.sum(adata_v.X[i*adata.n_obs : (i+1)*adata.n_obs], axis=0))
    df = pd.DataFrame(sum_cells).T.round(2)
    df.columns = layers
    df.index = adata.var.index
    df_norm = df.div(df.sum(axis=1), axis=0).round(2)
    return df_norm

def df_heatmap_clustermap(df, layers, heat=False, cluster=True, row_colors=None, zoom=False, thres=None, figsize=(5,5)):
    '''heatmap, clustermap for df of gene x n_layer'''
    # TODO add figure title
    plt.figure(figsize=figsize)
    if heat==True:
        sns.heatmap(df)
    if cluster==True:
        sns.clustermap(df, row_colors=row_colors, cmap='mako')
    if zoom==True:
        subset = []
        for layer in layers:
            subset.append(df.loc[df[layer]>thres])
        if thres <= 0.5:
            print('No.')
            return
        df_subset = pd.concat(subset, axis=0) #
        sns.heatmap(df_subset)
        sns.clustermap(df_subset)
        return df_subset
    return

def df_to_adata(df, obs_label=None):
    '''table contents to X; index to obs; column to var;
        obs_label -- array-like labels for obs, len == len(obs)'''
    obs = pd.DataFrame(df.index)
    if obs_label:
        obs['label'] = obs_label
        adata = sc.AnnData(X=df.to_numpy(), var=pd.DataFrame(df.columns), obs=obs)
    else: adata = sc.AnnData(X=df.to_numpy(), var=pd.DataFrame(df.columns), obs=obs)
    return adata

def dict_to_long(Dict, key_type, var_list, var_type, normalize=True):
    '''convert dictionary of dataframes into long-form'''
    Dict1 = copy(Dict)
    for key in Dict1.keys():
        if normalize:
            Dict1[key] = Dict1[key].div(Dict1[key].sum(axis=1), axis=0).round(2)
        Dict1[key][key_type] = [key] * Dict1[key].shape[0]
    df_all = pd.concat(Dict1.values())
    df_var = dict()
    for var in var_list:
        df_var[var] = df_all[[var, key_type]].rename(columns={var: 'value'})
        df_var[var][var_type] = [var] * df_var[var].shape[0]
    return pd.concat(df_var[var] for var in var_list)

def df_wide_to_long(df, cols, diff, val_name):
    '''formatted for seaborn --ughhhhh'''
    df_col = dict()
    for c in cols:
        df_col[c] = pd.DataFrame(df[c])
        df_col[c][diff] = c
        df_col[c].rename(columns={c:val_name}, inplace=True)
    df_long = pd.concat([df_col[c] for c in cols])
    return df_long

def cell_iter(rp, curr_cell, nuclei, reads_plot, div=[0,10]):
    '''pipeline-specific (not generalized yet). Input regionprops & nuclei mask (of entire image)'''
    # mask and bbox of current cell & nuclei
    cb = rp[curr_cell].bbox
    cell_mask = rp[curr_cell].image
    nuc_mask = nuclei[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]
    nuc_mask[nuc_mask != curr_cell+1] = 0

    # dt
    cdt = dt(cell_mask)
    ndt = dt(np.logical_not(nuc_mask))
    area = np.logical_xor(cell_mask, nuc_mask)
    rdt = np.zeros(cdt.shape)
    rdt[area] = ndt[area] / cdt[area]

    # region
    regions = [nuc_mask] # list
    for d in range(1,len(div)):
        regions.append(np.logical_and(rdt > div[d-1], rdt <= div[d]))
    regions.append(rdt > div[-1])

    # count reads in each region
    count = list()
    reads_bbox = reads_plot[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]
    for r in regions:
        count.append(np.sum(np.logical_and(r, reads_bbox)))

    return count

def cell_iter_by_gene(rp, curr_cell, nuclei, gene_list, reads_plots_by_gene, div=[0,10]):
    # NOTE see plot_reads_by_gene
    '''count each gene in each region of each cell
        return a list of lists (n_region x n_genes)'''
    # mask and bbox of current cell & nuclei
    cb = rp[curr_cell].bbox
    cell_mask = rp[curr_cell].image
    nuc_mask = nuclei[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]  # TODO import nuclei
    nuc_mask[nuc_mask != curr_cell+1] = 0

    # dt
    cdt = dt(cell_mask)
    ndt = dt(np.logical_not(nuc_mask))
    area = np.logical_xor(cell_mask, nuc_mask)
    rdt = np.zeros(cdt.shape)
    rdt[area] = ndt[area] / cdt[area]

    # region
    regions = [nuc_mask] # list
    for d in range(1,len(div)):
        regions.append(np.logical_and(rdt > div[d-1], rdt <= div[d]))
    regions.append(rdt > div[-1])

    # count reads in each region
    count = []
    for r in regions:
        r_count = []
        for gene in gene_list:
            reads_plot = reads_plots_by_gene[gene] # reads_plot of the current gene
            reads_bbox = reads_plot[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]
            r_count.append(np.sum(np.logical_and(r, reads_bbox)))
        count.append(r_count)
    return count

def plot_reads(points, bases, image_shape, genes_to_exclude=None, genes_to_include=None):
    '''plot reads on image from points (.mat)
        give a list of either include or exclude'''
    if genes_to_exclude:
        points_selected = points[np.isin(bases, genes_to_exclude, invert=True)]
    elif genes_to_include:
        points_selected = points[np.isin(bases, genes_to_include, invert=False)]
    else: points_selected = points

    plot = np.zeros(image_shape, dtype=np.bool_)
    plot[tuple(points_selected.T)] = True
    # for pt in points_selected:
    #     plot[pt[0], pt[1], pt[2]] = 1
    return plot

def plot_reads_by_gene(points, bases, image_shape, gene_list):
    # NOTE can cause memory problem if single image is too large
    '''return a dictionary of plots with genes as keys and value of plots, each plotting only one gene'''
    plots = dict()
    for gene in tqdm(gene_list):
        pts = points[np.isin(bases, gene)]
        plot = np.zeros(image_shape, dtype=np.bool_)
        plot[tuple(pts.T)] = True
        plots[gene] = plot
    return plots

def check_point_in_bbox_3D(point, region_bbox):
    '''useless for now'''
    if len(point) != 3 or len(region_bbox) != 6:
        print('Shape error.')
        return False
    if point[0] >= region_bbox[0] and point[0] < region_bbox[3]: # regionprops - half interval
        if point[1] >= region_bbox[1] and point[1] < region_bbox[4]:
            if point[2] >= region_bbox[2] and point[2] < region_bbox[5]:
                return True
    return False

def genes_to_index(genes):
    gti = dict()
    for i, gene in enumerate(genes):
        gti[gene] = i
    return gti

def cell_iter_new(rp, curr_cell, nuclei, points, bases, genes_to_index, div=[0,10]):
    ## TODO fix bugs
    '''multiple genes included; return a list of counts of length n_regions x 1 x n_genes as rows in layers of AnnData'''
    # mask and bbox of current cell & nuclei
    cb = rp[curr_cell].bbox
    cell_mask = rp[curr_cell].image
    nuc_mask = nuclei[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]
    nuc_mask[nuc_mask != curr_cell+1] = 0

    # dt
    cdt = dt(cell_mask)
    ndt = dt(np.logical_not(nuc_mask))
    area = np.logical_xor(cell_mask, nuc_mask)
    rdt = np.zeros(cdt.shape)
    rdt[area] = ndt[area] / cdt[area]

    # region
    regions = [nuc_mask] # list NOTE alternatively using label image
    for d in range(1,len(div)):
        regions.append(np.logical_and(rdt > div[d-1], rdt <= div[d]))
    regions.append(rdt > div[-1])

    # get reads and bases in the bbox
    cell_pts = []
    cell_bs = []
    for i in range(len(points)):
        if check_point_in_bbox_3D(points[i], cb): ## ERROR!!!
            cell_pts.append(points[i])
            cell_bs.append(bases[i]) # cell_pts, cell_bs still aligned

    # count reads in each region
    count = []
    for r in regions:
        r_count = [0] * len(genes_to_index) # 1 x n_genes
        for i in range(len(cell_pts)):
            pts = cell_pts[i]
            if r[pts[0]-cb[0], pts[1]-cb[1], pts[2]-cb[2]]:
                r_count[genes_to_index[cell_bs[i]]] += 1
    count.append(r_count) # n_region x 1 x n_genes

    return count

def cell_iter_continuous(rp, curr_cell, cell_seg, nuclei_seg, curr_points, sampling=[1,1,1]):
    '''return a list of distance-ratio for all cytosol reads of a particular genes in all cells'''
    # NOTE -- DO NOT use when running multiple genes!!! bc u do not want to run DT over all cells again

    # mask and bbox of current cell & nuclei
    cb = rp[curr_cell].bbox
    cell_mask = rp[curr_cell].image
    nuc_mask = nuclei_seg[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]
    nuc_mask[nuc_mask != curr_cell+1] = 0

    # dt
    cdt = dt(cell_mask, sampling=sampling)
    ndt = dt(np.logical_not(nuc_mask), sampling=sampling)
    area = np.logical_xor(cell_mask, nuc_mask)
    rdt = np.zeros(cdt.shape)
    rdt[area] = ndt[area] / cdt[area]

    # get reads in the current cell's cytosol
    cell_pts = [pt for pt in curr_points if ((cell_seg[pt[0], pt[1], pt[2]] == curr_cell+1) and (nuclei_seg[pt[0], pt[1], pt[2]] == 0))]

    # record DR for each read
    dr = []
    for pt in cell_pts:
        r = rdt[pt[0]-cb[0], pt[1]-cb[1], pt[2]-cb[2]]
        dr.append(r/(r+1))

    return dr

def cell_iter_continuous_ALL_GENES(rp, curr_cell, nuclei_seg, points, bases, reads_assign_cell, reads_assign_nucleus, collector, sampling=[1,1,1]):
    '''return a list of distance-ratio for all cytosol reads of a particular genes in all cells
        Collector -- provided to write result to, with struct: dict{keys: gene_list; vals: [...]}'''

    # mask and bbox of current cell & nuclei
    cb = rp[curr_cell].bbox
    cell_mask = rp[curr_cell].image
    nuc_mask = nuclei_seg[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]
    nuc_mask[nuc_mask != curr_cell+1] = 0

    # dt
    cdt = dt(cell_mask, sampling=sampling)
    ndt = dt(np.logical_not(nuc_mask), sampling=sampling)
    area = np.logical_xor(cell_mask, nuc_mask)
    rdt = np.zeros(cdt.shape)
    rdt[area] = ndt[area] / cdt[area]

    # get cytosol points & corr. bases of current cell
    cell_pts = points[np.logical_and(reads_assign_cell==curr_cell+1, reads_assign_nucleus==0)]
    cell_bs = bases[np.logical_and(reads_assign_cell==curr_cell+1, reads_assign_nucleus==0)]
    for i, bs in enumerate(cell_bs):
        pt = cell_pts[i]
        r = rdt[pt[0]-cb[0], pt[1]-cb[1], pt[2]-cb[2]]
        collector[bs].append(r/(r+1))
    return

def cell_iter_cont_all(rp, curr_cell, cell_index, nuclei_seg, points, bases, reads_assign_cell, reads_assign_nucleus, cell_by_barcode, genesToIndex, sampling=[1,1,1]):
    '''record into cell_by_barcode (collector), each entry (cell_x_gene) has a list of DR values'''
    # mask and bbox of current cell & nuclei
    cb = rp[curr_cell].bbox
    cell_mask = rp[curr_cell].image
    nuc_mask = nuclei_seg[cb[0]:cb[3], cb[1]:cb[4], cb[2]:cb[5]]
    nuc_mask[nuc_mask != curr_cell+1] = 0

    # dt
    cdt = dt(cell_mask, sampling=sampling)
    ndt = dt(np.logical_not(nuc_mask), sampling=sampling)
    area = np.logical_xor(cell_mask, nuc_mask)
    rdt = np.zeros(cdt.shape)
    rdt[area] = ndt[area] / cdt[area] # TODO see if where= solves

    # get cytosol points & corr. bases of current cell
    cell_pts = points[np.logical_and(reads_assign_cell==curr_cell+1, reads_assign_nucleus==0)]
    cell_bs = bases[np.logical_and(reads_assign_cell==curr_cell+1, reads_assign_nucleus==0)]
    for i, bs in enumerate(cell_bs):
        pt = cell_pts[i]
        r = rdt[pt[0]-cb[0], pt[1]-cb[1], pt[2]-cb[2]]
        # if type(cell_by_barcode[cell_index, genesToIndex[bs]]) == int: ## when initializing np.zeros(dtype=object) gives int (idk why but sounds unstable)
        #     cell_by_barcode[cell_index, genesToIndex[bs]] = [r/(r+1)]
        # else: cell_by_barcode[cell_index, genesToIndex[bs]].append(r/(r+1))
        if type(cell_by_barcode[cell_index, genesToIndex[bs]]) == list:
            cell_by_barcode[cell_index, genesToIndex[bs]].append(r/(r+1))
        else: cell_by_barcode[cell_index, genesToIndex[bs]] = [r/(r+1)]
    return

def compare_DR_distribution(sample_list, points_list, bases_list, rp_list, cell_list, nucleus_list, gene_list, sample=[1,1,1]):
    '''comparing DR distribution across samples'''
    # store as dictionaries with layer1keys - genes; layer2keys - samples
    dr = dict()
    for gene in gene_list:
        dr[gene] = dict()
        for sample in sample_list:
            dr[gene][sample] = []

    # process
    for gene in gene_list:
        for i, sample in enumerate(sample_list):
            curr_pts = points_list[i][bases_list[i]==gene]
            for cell in tqdm(range(len(rp_list[i]))):
                dr[gene][sample].extend(cell_iter_continuous(cell, cell_list[i], nucleus_list[i], curr_pts, rp=rp_list[i]))
    return dr

def adata_DRregion_as_layer(adata, cell, nuclei, bases, points, genes, div=[0,10], regions=['nuc', 'mid', 'peri'], save_adata=True, new_name='output_adata.h5ad'):
    '''add layers onto adata which represents gene expression in this subcellular region,
        defined by distance-ratio segmentation (partition).
        Returns a new AnnData with X, obs, var unchanged but layers as regions from compartments'''

    # rp, gti
    rp = measure.regionprops(cell)
    geneToIndex = genes_to_index(genes)

    # iter by cells
    cell_counts = []
    for curr_cell in tqdm(adata.obs['orig_index']):
        cell_counts.append(cell_iter_new(rp, curr_cell, nuclei, points, bases, geneToIndex, div))

    # create new AnnData
    layer = dict()
    for r in regions:
        layer[r] = []
    for ct in cell_counts:
        for i, re in enumerate(regions):
            layer[re].append(ct[i])
    new_adata = copy(adata)
    for r in layer:
        r_layer = np.array(layer[r])
        new_adata.layers[r] = layer[r]

    if save_adata:
        new_adata.write(new_name)
    return new_adata

def min_max(data):
    min_max_data = []
    x_min = min(data)
    x_max = max(data)
    for d in data:
        min_max_data.append((d-x_min)/(x_max-x_min))
    return min_max_data

def dr_norm(data):
    norm = []
    for d in data:
        norm.append(d/(d+1))
    return norm
