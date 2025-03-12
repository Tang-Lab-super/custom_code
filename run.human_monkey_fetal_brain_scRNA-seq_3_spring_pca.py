import os
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import sklearn.cluster
from sklearn.decomposition import PCA,TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csc_matrix
import scipy.io
import pickle
from scipy import sparse
from scipy.spatial.distance import pdist
from datetime import datetime
import json
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import matplotlib.pyplot as plt
import seaborn as sns

#########################################################################################################################################################
######################################################################### function ######################################################################
#########################################################################################################################################################
def tot_counts_norm(E, exclude_dominant_frac = 1, included = [], target_mean = 0):
    '''
    Cell-level total counts normalization of input counts matrix, excluding overly abundant genes if desired.
    Return normalized counts, average total counts, and (if exclude_dominant_frac < 1) list of genes used to calculate total counts
    '''

    E = E.tocsc()
    ncell = E.shape[0]
    if len(included) == 0:
        if exclude_dominant_frac == 1:
            tots_use = E.sum(axis=1)
        else:
            tots = E.sum(axis=1)
            wtmp = scipy.sparse.lil_matrix((ncell, ncell))
            wtmp.setdiag(1. / tots)
            included = np.asarray(~(((wtmp * E) > exclude_dominant_frac).sum(axis=0) > 0))[0,:]
            tots_use = E[:,included].sum(axis = 1)
            #print 'Excluded %i genes from normalization' %(np.sum(~included))
    else:
        tots_use = E[:,included].sum(axis = 1)

    if target_mean == 0:
        target_mean = np.mean(tots_use)

    w = scipy.sparse.lil_matrix((ncell, ncell))
    w.setdiag(float(target_mean) / tots_use)
    Enorm = w * E

    return Enorm.tocsc(), target_mean, included

def get_pca(E, base_ix=[], numpc=50, keep_sparse=False, normalize=True):
    '''
    Run PCA on the counts matrix E, gene-level normalizing if desired
    Return PCA coordinates
    '''
    # If keep_sparse is True, gene-level normalization maintains sparsity
    #     (no centering) and TruncatedSVD is used instead of normal PCA.

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    if keep_sparse:
        if normalize:
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply(E.T, 1 / zstd).T
        else:
            Z = E
        pca = TruncatedSVD(n_components=numpc)

    else:
        if normalize:
            zmean = E[base_ix,:].mean(0)
            zstd = np.sqrt(sparse_var(E[base_ix,:]))
            Z = sparse_multiply((E - zmean).T, 1/zstd).T
        else:
            Z = E
        pca = PCA(n_components=numpc)

    pca.fit(Z[base_ix,:])
    return pca.transform(Z)


def filter_genes(E, base_ix = [], min_vscore_pctl = 85, min_counts = 3, min_cells = 3, show_vscore_plot = False, sample_name = ''):
    '''
    Filter genes by expression level and variability
    Return list of filtered gene indices
    '''

    if len(base_ix) == 0:
        base_ix = np.arange(E.shape[0])

    Vscores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b = get_vscores(E[base_ix, :])
    ix2 = Vscores>0
    Vscores = Vscores[ix2]
    gene_ix = gene_ix[ix2]
    mu_gene = mu_gene[ix2]
    FF_gene = FF_gene[ix2]
    min_vscore = np.percentile(Vscores, min_vscore_pctl)
    ix = (((E[:,gene_ix] >= min_counts).sum(0).A.squeeze() >= min_cells) & (Vscores >= min_vscore))

    if show_vscore_plot:
        import matplotlib.pyplot as plt
        x_min = 0.5*np.min(mu_gene)
        x_max = 2*np.max(mu_gene)
        xTh = x_min * np.exp(np.log(x_max/x_min)*np.linspace(0,1,100))
        yTh = (1 + a)*(1+b) + b * xTh
        plt.figure(figsize=(8, 6));
        plt.scatter(np.log10(mu_gene), np.log10(FF_gene), c = [.8,.8,.8], alpha = 0.3, edgecolors='');
        plt.scatter(np.log10(mu_gene)[ix], np.log10(FF_gene)[ix], c = [0,0,0], alpha = 0.3, edgecolors='');
        plt.plot(np.log10(xTh),np.log10(yTh));
        plt.title(sample_name)
        plt.xlabel('log10(mean)');
        plt.ylabel('log10(Fano factor)');
        plt.show()

    return gene_ix[ix]

def get_vscores(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1):
    '''
    Calculate v-score (above-Poisson noise statistic) for genes in the input counts matrix
    Return v-scores and other stats
    '''

    ncell = E.shape[0]

    mu_gene = E.mean(axis=0).A.squeeze()
    gene_ix = np.nonzero(mu_gene > min_mean)[0]
    mu_gene = mu_gene[gene_ix]

    tmp = E[:,gene_ix]
    tmp.data **= 2
    var_gene = tmp.mean(axis=0).A.squeeze() - mu_gene ** 2
    del tmp
    FF_gene = var_gene / mu_gene

    data_x = np.log(mu_gene)
    data_y = np.log(FF_gene / mu_gene)

    x, y = runningquantile(data_x, data_y, fit_percentile, nBins)
    x = x[~np.isnan(y)]
    y = y[~np.isnan(y)]

    gLog = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])
    h,b = np.histogram(np.log(FF_gene[mu_gene>0]), bins=200)
    b = b[:-1] + np.diff(b)/2
    max_ix = np.argmax(h)
    c = np.max((np.exp(b[max_ix]), 1))
    errFun = lambda b2: np.sum(abs(gLog([x,c,b2])-y) ** error_wt)
    b0 = 0.1
    b = scipy.optimize.fmin(func = errFun, x0=[b0], disp=False)
    a = c / (1+b) - 1


    v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
    CV_eff = np.sqrt((1+a)*(1+b) - 1);
    CV_input = np.sqrt(b);

    return v_scores, CV_eff, CV_input, gene_ix, mu_gene, FF_gene, a, b

def runningquantile(x, y, p, nBins):

    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    dx = (x[-1] - x[0]) / nBins
    xOut = np.linspace(x[0]+dx/2, x[-1]-dx/2, nBins)

    yOut = np.zeros(xOut.shape)

    for i in range(len(xOut)):
        ind = np.nonzero((x >= xOut[i]-dx/2) & (x < xOut[i]+dx/2))[0]
        if len(ind) > 0:
            yOut[i] = np.percentile(y[ind], p)
        else:
            if i > 0:
                yOut[i] = yOut[i-1]
            else:
                yOut[i] = np.nan

    return xOut, yOut

def sparse_var(E, axis=0):
    ''' variance across the specified axis '''
    mean_gene = E.mean(axis=axis).A.squeeze()
    tmp = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene ** 2

def sparse_multiply(E, a):
    ''' multiply each row of E by a scalar '''
    nrow = E.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E

def sparse_zscore(E):
    ''' z-score normalize each column of E '''
    mean_gene = E.mean(0)
    stdev_gene = np.sqrt(sparse_var(E))
    return sparse_multiply((E - mean_gene).T, 1/stdev_gene).T

def get_knn_graph(X, k=5, dist_metric='euclidean', approx=False, return_edges=True):
    '''
    Build k-nearest-neighbor graph
    Return edge list and nearest neighbor matrix
    '''

    t0 = time.time()
    if approx:
        try:
            from annoy import AnnoyIndex
        except:
            approx = False
            #print 'Could not find library "annoy" for approx. nearest neighbor search'
    if approx:
        #print 'Using approximate nearest neighbor search'

        if dist_metric == 'cosine':
            dist_metric = 'angular'
        npc = X.shape[1]
        ncell = X.shape[0]
        annoy_index = AnnoyIndex(npc, metric=dist_metric)

        for i in range(ncell):
            annoy_index.add_item(i, list(X[i,:]))
        annoy_index.build(10) # 10 trees

        knn = []
        for iCell in range(ncell):
            knn.append(annoy_index.get_nns_by_item(iCell, k + 1)[1:])
        knn = np.array(knn, dtype=int)

    else:
        #print 'Using sklearn NearestNeighbors'

        if dist_metric == 'cosine':
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric, algorithm='brute').fit(X)
        else:
            nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)
        knn = nbrs.kneighbors(return_distance=False)

    if return_edges:
        links = set([])
        for i in range(knn.shape[0]):
            for j in knn[i,:]:
                links.add(tuple(sorted((i,j))))

        t_elapse = time.time() - t0
        #print 'kNN graph built in %.3f sec' %(t_elapse)

        return links, knn
    return knn

def write_color_tracks(ctracks, fname):
    out = []
    for name,score in ctracks.items():
        line = name + ',' + ','.join(['%.3f' %x for x in score])
        out += [line]
    out = sorted(out,key=lambda x: x.split(',')[0])
    open(fname,'w').write('\n'.join(out))

def get_color_stats_genes(color_stats, E, gene_list):
    means = E.mean(0).A.squeeze()
    stdevs = np.sqrt(sparse_var(E, 0))
    mins = E.min(0).todense().A1
    maxes = E.max(0).todense().A1

    pctl = 99.6
    pctl_n = (100-pctl) / 100. * E.shape[0]
    pctls = np.zeros(E.shape[1], dtype=float)
    for iG in range(E.shape[1]):
        n_nonzero = E.indptr[iG+1] - E.indptr[iG]
        if n_nonzero > pctl_n:
            pctls[iG] = np.percentile(E.data[E.indptr[iG]:E.indptr[iG+1]], 100 - 100 * pctl_n / n_nonzero)
        else:
            pctls[iG] = 0
        color_stats[gene_list[iG]] = tuple(map(float, (means[iG], stdevs[iG], mins[iG], maxes[iG], pctls[iG])))
    return color_stats

def get_color_stats_custom(color_stats, custom_colors):
    for k,v in custom_colors.items():
        color_stats[k] = tuple(map(float, (np.mean(v),np.std(v),np.min(v),np.max(v),np.percentile(v,99))))
    return color_stats

def save_color_stats(filename, color_stats):
    with open(filename,'w') as f:
        f.write(json.dumps(color_stats,indent=4, sort_keys=True))

def build_categ_colors(categorical_coloring_data, cell_groupings):
    for k,labels in cell_groupings.items():
        label_colors = {l:frac_to_hex(float(i)/len(set(labels))) for i,l in enumerate(list(set(labels)))}
        categorical_coloring_data[k] = {'label_colors':label_colors, 'label_list':labels}
    return categorical_coloring_data

def save_cell_groupings(filename, categorical_coloring_data):
    with open(filename,'w') as f:
        f.write(json.dumps(categorical_coloring_data,indent=4, sort_keys=True))

def write_graph(filename, n_nodes, edges):
    nodes = [{'name':int(i),'number':int(i)} for i in range(n_nodes)]
    edges = [{'source':int(i), 'target':int(j), 'distance':0} for i,j in edges]
    out = {'nodes':nodes,'links':edges}
    open(filename,'w').write(json.dumps(out,indent=4, separators=(',', ': ')))

def write_edges(filename, edges):
    with open(filename, 'w') as f:
        for e in edges:
            f.write('%i;%i\n' %(e[0], e[1]))

def save_spring_dir_sparse_hdf5(E,gene_list,project_directory, edges, custom_colors={}, cell_groupings={}):

    if not os.path.exists(project_directory):
        os.makedirs(project_directory)

    if not project_directory[-1] == '/':
        project_directory += '/'

    # save custom colors
    custom_colors['Uniform'] = np.zeros(E.shape[0])
    write_color_tracks(custom_colors, project_directory+'color_data_gene_sets.csv')

    # create and save a dictionary of color profiles to be used by the visualizer
    color_stats = {}
    color_stats = get_color_stats_genes(color_stats, E, gene_list)
    color_stats = get_color_stats_custom(color_stats, custom_colors)
    # return color_stats
    save_color_stats(project_directory + 'color_stats.json', color_stats)

    # save cell labels
    categorical_coloring_data = {}
    categorical_coloring_data = build_categ_colors(categorical_coloring_data, cell_groupings)
    save_cell_groupings(project_directory+'categorical_coloring_data.json', categorical_coloring_data)

    # write graph
    write_graph(project_directory + 'graph_data.json', E.shape[0], edges)
    write_edges(project_directory + 'edges.csv', edges)

def get_force_layout(links, n_cells, n_iter=100, edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, jitterTolerance=1, verbose=False):
    from fa2 import ForceAtlas2
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(n_cells))
    G.add_edges_from(list(links))

    forceatlas2 = ForceAtlas2(
                  # Behavior alternatives
                  outboundAttractionDistribution=False,  # Dissuade hubs
                  linLogMode=False,  # NOT IMPLEMENTED
                  adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                  edgeWeightInfluence=edgeWeightInfluence,

                  # Performance
                  jitterTolerance=jitterTolerance,  # Tolerance
                  barnesHutOptimize=True,
                  barnesHutTheta=barnesHutTheta,
                  multiThreaded=False,  # NOT IMPLEMENTED

                  # Tuning
                  scalingRatio=scalingRatio,
                  strongGravityMode=False,
                  gravity=gravity,
                  # Log
                  verbose=verbose)

    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=n_iter)
    positions = np.array([positions[i] for i in sorted(positions.keys())])
    return positions



save_path = '/data3/xiongdan/scRNA-seq/reRmDoublet/specificCells/d14_60_mergeProgenitor/'
os.chdir(save_path)
prefix = '5.harmony_embrddings'

dist_metric = 'euclidean'
use_approxnn = False

df_harmony = pd.read_csv('5.harmony_embrddings.txt', index_col=0, sep='\t')
df_anno = pd.read_csv('5.harmony_metadata_anno.txt', index_col=0, sep='\t')
# df_anno = df_anno[df_anno['SubAnnotation'].isin(['Plasma', 'B_Naive', 'B_Memory', 'B_GC', 'B_Other'])]
df_anno = df_anno[df_anno['SubAnnotation'].isin(['G2M', 'NEC', 'prog', 'IPC', 'EN', 'IN'])]
df_anno = df_anno.reindex(df_harmony.index)

Epca = df_harmony.values
# pipeline knn graph and links
# for k_neigh in range(30, 50):

for k_neigh in range(15,20,1):
    links_pipeline, knn_graph_pipeline = get_knn_graph(Epca, k=k_neigh, dist_metric = dist_metric, approx=use_approxnn)

    # custome knn graph
    dic_ct = dict(zip(df_anno.index, df_anno['SubAnnotation']))
    dic_ct_neighbor = {'G2M': ['G2M', 'NEC'], 'NEC': ['G2M', 'NEC', 'prog'], 'prog': ['NEC', 'prog', 'IPC'],
                       'IPC': ['prog', 'IPC', 'EN', 'IN'], 'EN': ['IPC', 'EN'], 'IN': ['IPC', 'IN']}

    nbrs = NearestNeighbors(n_neighbors=min(df_harmony.shape[0]-1, 5000), metric=dist_metric).fit(Epca)
    knn = nbrs.kneighbors(return_distance=False)
    # select cell type
    new_knn = []
    for idx_target, infor in enumerate(knn):
        new_neighbor = []
        target_celltype = dic_ct[df_harmony.index[idx_target]]
        right_neighbor_celltypes = dic_ct_neighbor[target_celltype]
        for idx_neighbor in infor:
            if len(new_neighbor) == k_neigh:
                break
            neighbor_celltype = dic_ct[df_harmony.index[idx_neighbor]]
            if neighbor_celltype in right_neighbor_celltypes:
                new_neighbor.append(idx_neighbor)
        new_knn.append(new_neighbor)

    knn_graph = np.vstack(new_knn)
    links = set([])
    for i in range(knn_graph.shape[0]):
        for j in knn_graph[i,:]:
            links.add(tuple(sorted((i,j))))


    num_force_iter = 1000
    if num_force_iter > 0:
        positions = get_force_layout(links, Epca.shape[0], n_iter=num_force_iter,
            edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05,
            jitterTolerance=1, verbose=False)
        positions = positions / 5.0
        positions = positions - np.min(positions, axis = 0) - np.ptp(positions, axis = 0) / 2.0
        positions[:,0] = positions[:,0]  + 750
        positions[:,1] = positions[:,1]  + 250

    df_position = np.hstack((np.arange(positions.shape[0])[:,None], positions))
    df_coor = pd.DataFrame(df_position[:, 1:3], index=df_anno.index, columns=['x', 'y'])


    dic_anno = dict(zip(df_anno.index, df_anno['SubAnnotation']))
    df_coor.index = df_harmony.index
    df_coor['SubAnnotation'] = df_anno.index.map(dic_anno)


    dic_color = {'G2M': '#5E4FA2', 'NEC': '#54AEAD', 'prog': '#BEE4A0', 'IPC':'#FDBE6F', 'EN':'#E95D46', 'IN':'#9E0142'}
    df_coor['color'] = df_coor['SubAnnotation'].map(dic_color)
    plt.figure()
    plt.scatter(x=df_coor['x'], y=df_coor['y'], c=df_coor['color'], s=0.5)
    plt.savefig(f'{save_path}/{prefix}.test.scatter_k{k_neigh}.pdf')


    import random
    def random_color():
        color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])]
        return color

    df_coor['sample'] = [x.split('_')[0] for x in df_coor.index]
    dic_color = {x:random_color()[0] for x in df_coor['sample'].unique()}
    df_coor['sample_color'] = df_coor['sample'].map(dic_color)
    plt.figure()
    plt.scatter(x=df_coor['x'], y=df_coor['y'], c=df_coor['sample_color'], s=0.5)
    plt.savefig(f'{save_path}/{prefix}.test.scatter.sample_k{k_neigh}.pdf')

    df_coor.to_csv(f'{save_path}/{prefix}.df_coor_k{k_neigh}.txt', sep='\t')

