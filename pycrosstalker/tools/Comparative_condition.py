import pandas as pd
import networkx as nx

def create_diff_table(data, out_path, comparison=None):
    """
    Read the lrobject and generate the comparative tables

    Parameters
    ----------
    data :
        LRObj with single condition
    out_path :
        output path
    
    Returns
    -------
    LRObject
    
    """

    def process_pair(exp_table, ctr_table):
        tmp_data = pd.merge(exp_table, ctr_table, on='allpair', how='outer')
        tmp_data[['ligpair', 'recpair']] = tmp_data['allpair'].str.split('@', expand=True)
        tmp_data['LRScore_x'] = tmp_data['LRScore_x'].fillna(0)
        tmp_data['LRScore_y'] = tmp_data['LRScore_y'].fillna(0)
        tmp_data['LRScore'] = tmp_data['LRScore_x'] - tmp_data['LRScore_y']
        
        final_data = tmp_data[tmp_data['LRScore'] != 0].copy()
        final_data['type_gene_A'] = final_data['type_gene_A_x'].combine_first(final_data['type_gene_A_y'])
        final_data['type_gene_B'] = final_data['type_gene_B_x'].combine_first(final_data['type_gene_B_y'])
        final_data['gene_A'] = final_data['gene_A_x'].combine_first(final_data['gene_A_y'])
        final_data['gene_B'] = final_data['gene_B_x'].combine_first(final_data['gene_B_y'])
        final_data['source'] = final_data['source_x'].combine_first(final_data['source_y'])
        final_data['target'] = final_data['target_x'].combine_first(final_data['target_y'])
        final_data['cellpair'] = final_data['source'] + "@" + final_data['target']
        final_data['interaction_type'] = final_data['type_gene_A'] + final_data['type_gene_B']
        
        final_data['interaction_type'] = final_data['interaction_type'].replace({
            'LigandReceptor': 'LR',
            'ReceptorTranscription Factor': 'RTF',
            'Transcription FactorLigand': 'TFL'
        })
        
        final = (final_data[~final_data['type_gene_A'].str.contains('Transcription') & 
                            ~final_data['type_gene_B'].str.contains('Transcription')]
                 .groupby('cellpair')['LRScore'].sum().reset_index())
        
        final[['u', 'v']] = final['cellpair'].str.split('@', expand=True)
        final.dropna(subset=['u', 'v'], inplace=True)
        
        raw_inter = final_data.loc[~final_data['type_gene_A'].str.contains('Transcription') & 
                                   ~final_data['type_gene_B'].str.contains('Transcription'),
                                   'cellpair'].value_counts()
        
        raw_inter = raw_inter.reindex(final['cellpair']).fillna(0)
        freq = (raw_inter - raw_inter.min()) / (raw_inter.max() - raw_inter.min()) + 0.1
        final['freq'] = freq.values
        final['pair'] = final['cellpair']
        
        return final_data, final

    if comparison is not None:
        for pair in comparison:
            ctr_name, exp_name = pair[1], pair[0]
            cmp_name = f"{exp_name}_x_{ctr_name}"
            exp_table, ctr_table = data['tables'][exp_name], data['tables'][ctr_name]
            final_data, final = process_pair(exp_table, ctr_table)
            data['tables'][cmp_name] = final_data
            
            G = nx.DiGraph()
            for _, row in final.iterrows():
                G.add_edge(row['u'], row['v'], LRScore=row['LRScore'], freq=row['freq'], weight=row['LRScore'], inter=row['freq'])
            
            data['graphs'][cmp_name] = G
            
            G_ggi = nx.DiGraph()
            for _, row in final_data.iterrows():
                G_ggi.add_edge(row['ligpair'], row['recpair'], LRScore=row['LRScore'], weight=row['LRScore'], inter=row['LRScore'])
            
            data['graphs_ggi'][cmp_name] = G_ggi
    else:
        ctr_name = list(data['tables'].keys())[0]
        ctr_table = data['tables'][ctr_name]
        for exp_name in list(data['tables'].keys())[1:]:
            cmp_name = f"{exp_name}_x_{ctr_name}"
            exp_table = data['tables'][exp_name]
            final_data, final = process_pair(exp_table, ctr_table)
            data['tables'][cmp_name] = final_data
            
            G = nx.DiGraph()
            for _, row in final.iterrows():
                G.add_edge(row['u'], row['v'], LRScore=row['LRScore'], freq=row['freq'], weight=row['LRScore'], inter=row['freq'])
            
            data['graphs'][cmp_name] = G
            
            G_ggi = nx.DiGraph()
            for _, row in final_data.iterrows():
                G_ggi.add_edge(row['ligpair'], row['recpair'], LRScore=row['LRScore'], weight=row['LRScore'], inter=row['LRScore'])
            
            data['graphs_ggi'][cmp_name] = G_ggi
    
    return data
