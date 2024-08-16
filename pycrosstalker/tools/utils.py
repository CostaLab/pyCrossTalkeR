import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import igraph
import itertools
from scipy.stats import fisher_exact
import pickle
from itertools import combinations
import os
import networkx as nx
import re

#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param data lrobject
#'@param out_path to save the lrobject with ranking
#'@param sel_columns columns to consider
#'@param slot slot of the networks graphs_ggi to gene cell interaction and abs
#'@import tibble
#'@import utils
#'@import dplyr
#'@return list
#'@importFrom tidyr %>%
#'@importFrom stats prcomp
#'@NoRd

def ranking(data, out_path, sel_columns, slot="graphs"):

    sc = StandardScaler()
    slot_data = data.get(slot, {})
    rankings = data.get('rankings', {})
    pca_data = data.get('pca', {})
    stats_data = data.get('stats', {})
    model = PCA(n_components= 2)
    for graph_name, graph in slot_data.items():
        if "_x_" in graph_name:  # Signed Analysis
            components = list(nx.connected_components(graph.to_undirected()))
            all_both = None
            
            for comp in components:
                subgraph = graph.subgraph(comp)
                if all_both is None:
                    all_both = ranking_net(subgraph, mode=False)
                else:
                    tmp = ranking_net(subgraph, mode=False)
                    all_both = pd.concat([all_both, tmp], ignore_index=True)
            
            if "_ggi" in slot:
                all_both = comparative_pagerank(rankings, slot, graph_name, all_both)
                all_both = comparative_med(rankings, slot, graph_name, all_both)
                rankings[f"{graph_name}_ggi"] = all_both
                pca_data_frame = all_both.drop(columns=['nodes']).loc[:, (all_both != 0).any(axis=0)]
                pca_data_frame = sc.fit_transform(pca_data_frame)
                pca = model.fit_transform(pca_data_frame)
                pca_data[f"{graph_name}_ggi"] = pca
            else:
                all_both = comparative_pagerank(rankings, slot, graph_name, all_both)
                all_both = comparative_med(rankings, slot, graph_name, all_both)
                rankings[graph_name] = all_both
                pca_data_frame = all_both.drop(columns=['nodes']).loc[:, (all_both != 0).any(axis=0)]
                pca_data_frame = sc.fit_transform(pca_data_frame)
                pca = model.fit_transform(pca_data_frame)
                pca_data[graph_name] = pca
        else:  # Unsigned Analysis
            components = list(nx.connected_components(graph.to_undirected()))
            all_data = None
            
            for comp in components:
                subgraph = graph.subgraph(comp)
                if all_data is None:
                    all_data = ranking_net(subgraph)
                else:
                    tmp = ranking_net(subgraph)
                    all_data = pd.concat([all_data, tmp], ignore_index=True)
            
            if "_ggi" in slot:
                final_data = None
                table = data['tables'][graph_name]
                cls = pd.unique(pd.concat([table[sel_columns[0]], table[sel_columns[1]]]))
                for i in cls:
                    all_eq = pd.unique(pd.concat([
                        table['ligpair'][table[sel_columns[0]] == i],
                        table['recpair'][table[sel_columns[1]] == i]
                    ]))
                    edges = pd.DataFrame(list(itertools.combinations(all_eq, 2)), columns=['u', 'v'])
                    edges['LRScore'] = 0.0
                    if final_data is None:
                        final_data = edges
                    else:
                        final_data = pd.concat([final_data, edges], ignore_index=True)
                
                tmp_tbl = table[['ligpair', 'recpair', 'LRScore']]
                all1 = pd.concat([tmp_tbl, final_data])
                tmp_net = nx.from_pandas_edgelist(all1, source='ligpair', target='recpair', edge_attr='LRScore', create_using=nx.DiGraph())
                pg = nx.pagerank(tmp_net, weight='LRScore')
                all_data['Pagerank'] = all_data['nodes'].map(pg)
                rankings[f"{graph_name}_ggi"] = all_data
                pca_data_frame = all_data.drop(columns=['nodes']).loc[:, (all_data != 0).any(axis=0)]
                pca = model.fit_transform(pca_data_frame)
                pca_data[f"{graph_name}_ggi"] = pca 
            else:
                pg = nx.pagerank(graph, weight='LRScore')
                all_data['Pagerank'] = all_data['nodes'].map(pg)
                rankings[graph_name] = all_data
                pca_data_frame = all_data.drop(columns=['nodes']).loc[:, (all_data != 0).any(axis=0)]
                pca_data_frame = sc.fit_transform(pca_data_frame)
                pca = model.fit_transform(pca_data_frame)
                pca_data[graph_name] = pca

    data['rankings'] = rankings
    data['pca'] = pca_data
    data['stats'] = stats_data

    with open(os.path.join(out_path, "LR_data_final.pkl"), "wb") as f:
        pickle.dump(data, f)
    
    return data



#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param rankings tables lrobject
#'@param slotname slot of the networks graphs_ggi to gene cell interaction and abs
#'@param graphname graph comparison name
#'@param curr.rkg ranking table
#'@importFrom tidyr %>%
#'@import stringr
#'@return list
#'@NoRd
#' Network Ranking method
#'
#'@param graph lrobject
#'@param mode  is TRUE if is comparive mode
#'@return list
#'@import igraph
#'@importFrom tidyr %>%
#'@NoRd

def ranking_net(graph, mode=True):
    nodes = list(graph.nodes)
    
    if not mode:
        # Positive weights subgraph
        pos_edges = [(u, v) for u, v, d in graph.edges(data=True) if d['weight'] > 0]
        pos_subgraph = graph.edge_subgraph(pos_edges)
        deg_in_pos = dict(pos_subgraph.in_degree()) if isinstance(graph, nx.DiGraph) else dict(pos_subgraph.degree())
        deg_out_pos = dict(pos_subgraph.out_degree()) if isinstance(graph, nx.DiGraph) else dict(pos_subgraph.degree())
        
        # Negative weights subgraph
        neg_edges = [(u, v) for u, v, d in graph.edges(data=True) if d['weight'] < 0]
        neg_subgraph = graph.edge_subgraph(neg_edges)
        deg_in_neg = dict(neg_subgraph.in_degree()) if isinstance(graph, nx.DiGraph) else dict(neg_subgraph.degree())
        deg_out_neg = dict(neg_subgraph.out_degree()) if isinstance(graph, nx.DiGraph) else dict(neg_subgraph.degree())
        
        # Ensure nodes exist in the dictionaries
        deg_in_pos = {node: deg_in_pos.get(node, 0) + 1 for node in nodes}
        deg_out_pos = {node: deg_out_pos.get(node, 0) + 1 for node in nodes}
        deg_in_neg = {node: deg_in_neg.get(node, 0) + 1 for node in nodes}
        deg_out_neg = {node: deg_out_neg.get(node, 0) + 1 for node in nodes}
        
        centrality_table = pd.DataFrame({
            'nodes': nodes,
            'Listener': [round(deg_in_pos[node] - deg_in_neg[node], 2) for node in nodes],
            'Influencer': [round(deg_out_pos[node] - deg_out_neg[node], 2) for node in nodes]
        })
    else:
        # Calculating degrees
        deg_in_pos = dict(graph.in_degree()) if isinstance(graph, nx.DiGraph) else dict(graph.degree())
        deg_out_pos = dict(graph.out_degree()) if isinstance(graph, nx.DiGraph) else dict(graph.degree())
        
        # Calculating betweenness centrality with absolute weights
        abs_weights_graph = graph.copy()
        for u, v, d in abs_weights_graph.edges(data=True):
            d['weight'] = abs(d['weight'])
        bet = nx.betweenness_centrality(abs_weights_graph, weight='weight', normalized= False)
        
        centrality_table = pd.DataFrame({
            'nodes': nodes,
            'Listener': [round(deg_in_pos[node], 2) for node in nodes],
            'Influencer': [round(deg_out_pos[node], 2) for node in nodes],
            'Mediator': [round(bet[node], 2) for node in nodes]
        })
    
    centrality_table.fillna(0, inplace=True)
    return centrality_table


#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param rankings tables lrobject
#'@param slotname slot of the networks graphs_ggi to gene cell interaction and abs
#'@param graphname graph comparison name
#'@param curr.rkg ranking table
#'@importFrom tidyr %>%
#'@import stringr
#'@return list
#'@NoRd

def comparative_pagerank(rankings, slotname, graphname, curr_rkg):
    p_f1 = p_f2 = 0.5  # probability to be at disease
    allnodes = pd.DataFrame(curr_rkg['nodes'], columns=['nodes'])
    if '_filtered' in graphname:
        curr = re.sub('_filtered', '', graphname)
        curr = re.split('_x_', curr)
        p_ctr = curr[1]
        q_exp = curr[0]
    else:
        curr = re.split('_x_', graphname)
        p_ctr = curr[1]
        q_exp = curr[0]

 
    if "_ggi" in slotname:
        p = rankings[p_ctr + '_ggi'][['nodes','Pagerank']]
        q = rankings[q_exp + '_ggi'][['nodes','Pagerank']]

        

    else:
        p = rankings[p_ctr][['nodes', 'Pagerank']].loc[rankings[q_exp]['Pagerank'].index]
        q = rankings[q_exp][['nodes', 'Pagerank']]

    
    p.columns = ['nodes', 'p_ctr']
    q.columns = ['nodes', 'p_dis']

    final = pd.merge(allnodes,p, on='nodes', how='left')
    final = pd.merge(final, q, on='nodes', how='left')

    final['p_ctr'] = final['p_ctr'].fillna(0)
    final['p_dis'] = final['p_dis'].fillna(0)
    
    alpha = 0.01
    final['p_ctr'] = final['p_ctr'] + alpha
    final['p_dis'] = final['p_dis'] + alpha

    final['p_ctr'] = final['p_ctr'] / final['p_ctr'].sum()
    final['p_dis'] = final['p_dis'] / final['p_dis'].sum()

    p = final['p_ctr']
    q = final['p_dis']

    pc = p * p_f1 + q * p_f2
    pcontrol = (p_f1 * p) / pc
    pdisease = (p_f2 * q) / pc

    final_result = np.log(pdisease / pcontrol)  
    
    curr_rkg['Pagerank'] = final_result
    return curr_rkg


#'Delta betweenness the most interactive gene (ligand or receptor)
#'
#'@param rankings tables lrobject
#'@param slotname slot of the networks graphs_ggi to gene cell interaction and abs
#'@param graphname graph comparison name
#'@param curr.rkg ranking table
#'@importFrom tidyr %>%
#'@import stringr
#'@return list
#'@NoRd


def comparative_med(rankings, slotname, graphname, curr_rkg):
    allnodes = curr_rkg['nodes']
    curr = graphname.split('_x_')
    p_ctr = curr[1]
    q_exp = curr[0]

   
    if "_ggi" in slotname:
        p = rankings[p_ctr + '_ggi']['Mediator']
        q = rankings[q_exp + '_ggi']['Mediator']
    else:
        p = rankings[p_ctr]['Mediator'][rankings[q_exp]['Mediator'].index]
        q = rankings[q_exp]['Mediator']
    

    final = pd.DataFrame({'p.ctr': p, 'p.dis': q, 'names': allnodes})
    final['p.ctr'] = final['p.ctr'].fillna(0)
    final['p.dis'] = final['p.dis'].fillna(0)
    curr_rkg['Mediator'] = final['p.dis'] - final['p.ctr']
    return curr_rkg





def add_node_type(df):
    df['gene_A'] = df.apply(lambda row: f"{row['gene_A']}|L" if row['type_gene_A'] == "Ligand" else row['gene_A'], axis=1)
    df['gene_A'] = df.apply(lambda row: f"{row['gene_A']}|R" if row['type_gene_A'] == "Receptor" else row['gene_A'], axis=1)
    df['gene_A'] = df.apply(lambda row: f"{row['gene_A']}|TF" if row['type_gene_A'] == "Transcription Factor" else row['gene_A'], axis=1)
    
    df['gene_B'] = df.apply(lambda row: f"{row['gene_B']}|L" if row['type_gene_B'] == "Ligand" else row['gene_B'], axis=1)
    df['gene_B'] = df.apply(lambda row: f"{row['gene_B']}|R" if row['type_gene_B'] == "Receptor" else row['gene_B'], axis=1)
    df['gene_B'] = df.apply(lambda row: f"{row['gene_B']}|TF" if row['type_gene_B'] == "Transcription Factor" else row['gene_B'], axis=1)
    
    return df
