import pandas as pd
import numpy as np
import networkx as nx
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from anndata import AnnData
from .utils import *

def read_lr_single_condition(input, sel_columns, out_path="/tmp/", sep=",", colors=None):
    """
    This function loads the single conditions LR outputs and use it to generate the report data and it`s object It assumes that the table presents the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another measure

    Parameters
    ----------
    input :
        Named vector with the lrpaths of each output or an AnnData object
    sel_columns :
        selected columns
    out_path :
        Path to deposit the results
    sep :
        character used to divide the columns on input file
    colors :
        colorlist
    
    Returns
    -------
    LRObject
    
    """
    if isinstance(input, dict):
        annData = AnnData(uns = {'pycrosstalker': {'path': {}}})
        for cond, lrpath in input.items():
            # Reading data
            if isinstance(lrpath, str):
                data1 = pd.read_csv(lrpath, sep=sep)
                annData.uns['pycrosstalker']['path'][cond] = data1
            else:
                raise ValueError("Issue with input paths, please check!")
    elif isinstance(input, AnnData):
        annData = input
    else:
        raise ValueError("Input parameter must be either a file path or an AnnData object")

    data = {}
    graphs = {}
    graphs_ggi = {}
    unif_celltypes = []

    for cond in annData.uns['pycrosstalker']['path'].keys():
        data1 = annData.uns['pycrosstalker']['path'][cond]
        if not (data1['gene_A'].str.contains(r'\|').sum() > 0):
            data1 = add_node_type(data1)
        
        # Uniting columns
        data1['cellpair'] = data1[sel_columns[0]] + '@' + data1[sel_columns[1]]
        if 'Ligand' in data1[sel_columns[4]].unique():
            data1['ligpair'] = data1[sel_columns[0]] + '/' + data1[sel_columns[2]]
            data1['recpair'] = data1[sel_columns[1]] + '/' + data1[sel_columns[3]]
            data1['allpair'] = data1['ligpair'] + '@' + data1['recpair']
            data1[['ligpair', 'recpair']] = data1['allpair'].str.split('@', expand=True)
        else:
            data1[sel_columns[4]], data1[sel_columns[5]] = data1[sel_columns[5]], data1[sel_columns[4]]
            data1[sel_columns[3]], data1[sel_columns[4]] = data1[sel_columns[4]], data1[sel_columns[3]]
            data1['ligpair'] = data1[sel_columns[0]] + '/' + data1[sel_columns[2]]
            data1['recpair'] = data1[sel_columns[1]] + '/' + data1[sel_columns[3]]
            data1['allpair'] = data1['ligpair'] + '@' + data1['recpair']
            data1[['ligpair', 'recpair']] = data1['allpair'].str.split('@', expand=True)

        unif_celltypes += data1[sel_columns[0]].unique().tolist()
        unif_celltypes += data1[sel_columns[1]].unique().tolist()

        data1['LRScore'] = data1[sel_columns[-1]]

        final = data1.groupby('cellpair')['LRScore'].sum().reset_index()

        aux = final['cellpair'].str.split('@', expand=True)
        final['pair'] = final['cellpair']
        final[['u', 'v']] = aux

        filtervar = data1[sel_columns[4]].str.contains('Transcription') | data1[sel_columns[5]].str.contains('Transcription')
        raw_inter = data1.loc[~filtervar, 'cellpair'].value_counts()
        freq = (raw_inter[final['pair']].values - min(raw_inter)) / (max(raw_inter) - min(raw_inter)) + 0.1
        final['freq'] = freq

        # Convert to NetworkX
        graph1 = nx.DiGraph()
        for index, row in final.iterrows():
            graph1.add_edge(row['u'], 
                            row['v'], 
                            LRScore=row['LRScore'], 
                            freq=row['freq'], 
                            weight=row['LRScore'], 
                            inter=row['freq'])  # Add thickness (inter)

        graph2 = nx.DiGraph()
        for index, row in data1.iterrows():
            graph2.add_edge(row['ligpair'], 
                            row['recpair'], 
                            LRScore=row['LRScore'], 
                            mean=row['LRScore'], 
                            weight=row['LRScore'], 
                            inter=row['LRScore'])  # Add thickness (inter)

        data[cond] = data1
        graphs[cond] = nx.to_pandas_edgelist(graph1)
        graphs_ggi[cond] = nx.to_pandas_edgelist(graph2)

    # Create a full graph
    template = nx.complete_graph(len(set(unif_celltypes)))
    c  = nx.circular_layout(template)

    coords = {node: c[i] for i, node in enumerate(sorted(set(unif_celltypes), key=lambda x: x.lower()), start=0)}

    if colors is None:
        matplot_colors = plt.cm.get_cmap('Paired', len(set(unif_celltypes)))
        colors = {node: mcolors.to_hex(matplot_colors(i)) for i, node in enumerate(sorted(set(unif_celltypes), key=lambda x: x.lower()), start=0)}

    lr = {"graphs": graphs,
          "graphs_ggi": graphs_ggi,
          "tables": data,
          "colors": colors, 
          "coords": coords,
          "rankings": {},
          "pca" : {},
          "stats" : {}} 

    annData.uns['pycrosstalker']['results'] = lr
    return annData
