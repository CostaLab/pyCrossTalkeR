import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import igraph
import itertools
from scipy.stats import fisher_exact, MonteCarloMethod, mannwhitneyu
from scipy.cluster.hierarchy import linkage, leaves_list
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
    """
    Ranking the most interactive gene (ligand or receptor)

    Parameters
    ----------
    data :
        lrobject
    out_path :
        to save the lrobject with ranking
    sel_columns :
        columns to consider
    slot :
        slot of the networks graphs_ggi to gene cell interaction and abs
    
    Returns
    -------
    list
    
    """

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
    """
    Network Ranking method

    Parameters
    ----------
    graph :
        lrobject
    mode :
        is TRUE if is comparive mode
    
    Returns
    -------
    list
    
    """

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
    """
    Ranking the most interactive gene (ligand or receptor)

    Parameters
    ----------
    rankings :
        tables lrobject
    slotname :
        slot of the networks graphs_ggi to gene cell interaction and abs
    graphname :
        graph comparison name
    curr.rkg :
        ranking table
    
    Returns
    -------
    list
    
    """
    
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
    """
    Delta betweenness the most interactive gene (ligand or receptor)

    Parameters
    ----------
    rankings :
        tables lrobject
    slotname :
        slot of the networks graphs_ggi to gene cell interaction and abs
    graphname :
        graph comparison name
    curr.rkg :
        ranking table
    
    Returns
    -------
    list
    
    """
     
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
        p = rankings[p_ctr + '_ggi'][['nodes', 'Mediator']]
        q = rankings[q_exp + '_ggi'][['nodes', 'Mediator']]
    else:
        p = rankings[p_ctr][['nodes', 'Mediator']].loc[rankings[q_exp][['nodes', 'Mediator']].index]
        q = rankings[q_exp][['nodes', 'Mediator']]

    p.columns = ['nodes', 'm_ctr']
    q.columns = ['nodes', 'm_exp']

    final = pd.merge(allnodes,p, on='nodes', how='left')
    final = pd.merge(final, q, on='nodes', how='left')

    final['m_ctr'] = final['m_ctr'].fillna(0)
    final['m_exp'] = final['m_exp'].fillna(0)

    curr_rkg['Mediator'] = final['m_exp'] - final['m_ctr']
    return curr_rkg


def add_node_type(df):
    """
    Adding genetype to the gene names to distinguish biological function

    Parameters
    ----------
    df :
        dataframe with interaction data
    
    Returns
    -------
    df
    
    """

    df['gene_A'] = df.apply(lambda row: f"{row['gene_A']}|L" if row['type_gene_A'] == "Ligand" else row['gene_A'], axis=1)
    df['gene_A'] = df.apply(lambda row: f"{row['gene_A']}|R" if row['type_gene_A'] == "Receptor" else row['gene_A'], axis=1)
    df['gene_A'] = df.apply(lambda row: f"{row['gene_A']}|TF" if row['type_gene_A'] == "Transcription Factor" else row['gene_A'], axis=1)
    
    df['gene_B'] = df.apply(lambda row: f"{row['gene_B']}|L" if row['type_gene_B'] == "Ligand" else row['gene_B'], axis=1)
    df['gene_B'] = df.apply(lambda row: f"{row['gene_B']}|R" if row['type_gene_B'] == "Receptor" else row['gene_B'], axis=1)
    df['gene_B'] = df.apply(lambda row: f"{row['gene_B']}|TF" if row['type_gene_B'] == "Transcription Factor" else row['gene_B'], axis=1)
    
    return df


def fisher_test_cci(data, measure, out_path, comparison=None):
    """
    Evaluate Differences in the edge proportion

    Parameters
    ----------
    data :
        data
    measure :
        intensity
    out_path :
        save path
    
    Returns
    -------
    data (lr_object) with fisher stats
    
    """

    if comparison is not None:
        for pair in comparison:
            ctr_name = pair[1]
            exp_name = pair[0]
            
            c = data['tables'][ctr_name].groupby('cellpair').size().reset_index(name='measure')
            e = data['tables'][exp_name].groupby('cellpair').size().reset_index(name='measure')

            joined = pd.merge(c, e, on='cellpair', how='outer', suffixes=('_ctr', '_exp'))
            
            pvals = []
            measure_ctr_sum = joined['measure_ctr'].sum()
            measure_exp_sum = joined['measure_exp'].sum()
            for idx, row in joined.iterrows():
                ctotal = joined['measure_ctr'].sum() - row['measure_ctr']
                etotal = joined['measure_exp'].sum() - row['measure_exp']
                matrix = np.array([[row['measure_exp'], etotal], [row['measure_ctr'], ctotal]])

                n_resamples = 1000
                np.random.seed(42)  # For reproducibility

                matrix_data = np.hstack([
                    np.repeat(0, matrix[0, 0]),
                    np.repeat(1, matrix[0, 1]),
                    np.repeat(2, matrix[1, 0]),
                    np.repeat(3, matrix[1, 1]),
                ])

                bootstrap_p_values = []
                for _ in range(n_resamples):
                    resampled_data = np.random.choice(matrix_data, size=matrix_data.size, replace=True)
                    resampled_table = np.array([
                        [np.sum(resampled_data == 0), np.sum(resampled_data == 1)],
                        [np.sum(resampled_data == 2), np.sum(resampled_data == 3)]
                    ])
                    _, p_value = fisher_exact(resampled_table)
                    bootstrap_p_values.append(p_value)

                # Estimate the empirical p-value
                bootstrap_p_values = np.array(bootstrap_p_values)
                empirical_p_value = np.mean(bootstrap_p_values <= p_value)

                # Fisher's exact test
                _, p_value_1 = fisher_exact(matrix)

                # print(p_value_1)
                # print(empirical_p_value)
                
                # Log odds ratio
                odds_ratio = row['measure_exp'] / row['measure_ctr'] if row['measure_ctr'] != 0 else np.nan
                lodds = np.log2(odds_ratio) if odds_ratio > 0 else np.nan
                
                pvals.append({
                    'cellpair': row['cellpair'],
                    'p_value': p_value_1,
                    'lodds': lodds
                })

            pval_df = pd.DataFrame(pvals)
            data['stats'][f'{exp_name}_x_{ctr_name}'] = pval_df

        with open(os.path.join(out_path, "LR_data_final.pkl"), "wb") as f:
            pickle.dump(data, f)

        return data

    else:
        if len(data['tables']) >= 2:
            c_key = list(data['tables'].keys())[0]
            c = data['tables'][c_key].groupby('cellpair').size().reset_index(name='measure')

            for i in range(1, len(data['tables'])):
                if '_x_' not in list(data['tables'].keys())[i]:
                    e = data['tables'][list(data['tables'].keys())[i]].groupby('cellpair').size().reset_index(name='measure')

                    # Merge control and experimental data on 'cellpair'
                    joined = pd.merge(c, e, on='cellpair', how='inner', suffixes=('_ctr', '_exp'))
                    
                    pvals = []
                    measure_ctr_sum = joined['measure_ctr'].sum()
                    measure_exp_sum = joined['measure_exp'].sum()
                    for idx, row in joined.iterrows():
                        ctotal = measure_ctr_sum - row['measure_ctr']
                        etotal = measure_exp_sum - row['measure_exp']
                        matrix = np.array([[row['measure_exp'], etotal], [row['measure_ctr'], ctotal]])
                        odds_ratio, p_value = fisher_exact(matrix, alternative="two-sided")
    
                        # Monte Carlo resampling (if B > 0)
                        B = 1000
                        np.random.seed(42)  # For reproducibility
                        if B > 0:
                            count = 0
                            for _ in range(B):
                                shuffled = np.random.permutation(matrix.flatten()).reshape(matrix.shape)
                                if fisher_exact(shuffled, alternative="two-sided")[1] <= p_value:
                                    count += 1
                            perm_p_value = (count + 1) / (B + 1)  # Avoid zero probability
                        else:
                            perm_p_value = None

                        lodds = np.log2(odds_ratio) if odds_ratio > 0 else None  # Compute log odds ratio

                        pvals.append({
                            'cellpair': row['cellpair'],
                            'p_value': p_value,
                            'lodds': lodds
                        })

                    pval_df = pd.DataFrame(pvals)
                    data['stats'][f'{list(data["tables"].keys())[i]}_x_{list(data["tables"].keys())[0]}'] = pval_df

            with open(os.path.join(out_path, "LR_data_final.pkl"), "wb") as f:
                pickle.dump(data, f)

            return data


def filtered_graphs(data, out_path):
    """
    Filter comparison graphs by statistics

    Parameters
    ----------
    data :
        data
    out_path :
        save path
    
    Returns
    -------
    data (lr_object) with fisher comparison graphs
    
    """

    temp = {}
    for name, graph in data['graphs'].items():
        if '_x_' in name:
            h = [edge[0] for edge in graph.edges()]
            f = [edge[1] for edge in graph.edges()]

            stats = data['stats'][name]
            
            significant_edges = stats[stats['p_value'] <= 0.05]
            
            significant_edge_pairs = [
                (h[i], f[i]) for i, edge_pair in enumerate(zip(h, f))
                if f"{h[i]}@{f[i]}" in significant_edges['cellpair'].values
            ]

            filtered_graph = graph.edge_subgraph(significant_edge_pairs).copy()
            temp[name] = filtered_graph

    for name in temp:   
        data['graphs'][f"{name}_filtered"] = temp[name]

    with open(os.path.join(out_path, "LR_data_final.pkl"), "wb") as f:
        pickle.dump(data, f)

    return data


def mannwhitneyu_test_cci(data, measure, out_path, comparison=None):
    """
    Evaluate Differences in the edge strength

    Parameters
    ----------
    data :
        data
    measure :
        intensity
    out_path :
        save path
    
    Returns
    -------
    data (lr_object) with mannwittu stats
    
    """

    lcellpair = {}
    for key, df in data['tables'].items():
        lcellpair[key] = df['cellpair'].unique()

    if comparison:
        for pair in comparison:
            ctr_name, exp_name = pair

            results = []
            for cellpair in np.unique(np.concatenate(list(lcellpair.values()))):
                c = data['tables'][ctr_name].loc[data['tables'][ctr_name]['cellpair'] == cellpair, ['allpair', measure]]
                e = data['tables'][exp_name].loc[data['tables'][exp_name]['cellpair'] == cellpair, ['allpair', measure]]

                merged = pd.merge(c, e, on='allpair', how='outer').fillna(0)
                
                stat, p_value = mannwhitneyu(merged[measure + '_x'], merged[measure + '_y'], alternative='two-sided')
                
                results.append({
                    'cellpair': cellpair,
                    'statistic': stat,
                    'p_value': p_value
                })

            data['stats'][f'{exp_name}_x_{ctr_name}:MannU'] = pd.DataFrame(results)
    
    else:
        c_key = list(data['tables'].keys())[0]
        c = data['tables'][c_key]
        
        for key, df in data['tables'].items():
            if key != c_key and "_x_" not in key:

                results = []
                for cellpair in np.unique(np.concatenate(list(lcellpair.values()))):
                    c = data['tables'][c_key].loc[data['tables'][c_key]['cellpair'] == cellpair, ['allpair', measure]]
                    e = df.loc[df['cellpair'] == cellpair, ['allpair', measure]]

                    merged = pd.merge(c, e, on='allpair', how='outer').fillna(0)

                    stat, p_value = mannwhitneyu(merged[measure + '_x'], merged[measure + '_y'], alternative='two-sided')
                    eps = 1e-6
                    lfc = np.log2((merged[measure + '_y'].mean() + eps) / (merged[measure + '_x'].mean() + eps)) if merged[measure + '_x'].mean() > 0 else np.nan

                    results.append({
                        'cellpair': cellpair,
                        'statistic': stat,
                        'p_value': p_value,
                        'lfc': lfc
                    })

                data['stats'][f'{key}_x_{c_key}:MannU'] = pd.DataFrame(results)

    return data

def get_clustered_node_order(graph, weight="LRScore"):
    adj_df = nx.to_pandas_adjacency(graph, weight=weight).fillna(0)
    Z = linkage(adj_df.values, method="ward", metric="euclidean")
    row_order = leaves_list(Z)
    ordered_nodes = adj_df.index[row_order].tolist()
    return ordered_nodes

def create_ordered_circular_layout(ordered_nodes):
    n = len(ordered_nodes)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    layout = {
        node: (np.cos(angle), np.sin(angle))
        for node, angle in zip(ordered_nodes, angles)
    }
    return layout
