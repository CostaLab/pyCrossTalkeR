from pycrosstalker import tools as cttl
from pycrosstalker import plots as ctpl
import pickle
import pandas as pd

# Load data from pickle file
with open("../output/LR_data.pkl", "rb") as f:
    data = pickle.load(f)

print(f"\nTesting if data from pyCrossTalkeR is similar to data from CrossTalkeR")

def test_check_table_data():
    # Load table data from R
    data_R_CTR = pd.read_csv('R_data/table_CTR.csv')
    data_R_EXP = pd.read_csv('R_data/table_EXP.csv')
    data_R_EXP_x_CTR = pd.read_csv('R_data/table_EXP_x_CTR.csv')

    common_columns_CTR = data['tables']['CTR'].columns.intersection(data_R_CTR.columns)
    common_columns_EXP = data['tables']['EXP'].columns.intersection(data_R_EXP.columns)
    common_columns_EXP_x_CTR = data['tables']['EXP_x_CTR'].columns.intersection(data_R_EXP_x_CTR.columns)

    df1_sorted = data_R_EXP_x_CTR.sort_values(by="LRScore")[common_columns_EXP_x_CTR].reset_index(drop=True)
    df2_sorted = data['tables']['EXP_x_CTR'].sort_values(by="LRScore")[common_columns_EXP_x_CTR].reset_index(drop=True)

    df1_sorted['LRScore'] = df1_sorted['LRScore'].round(6)
    df2_sorted['LRScore'] = df2_sorted['LRScore'].round(6)

    assert data['tables']['CTR'][common_columns_CTR].equals(data_R_CTR[common_columns_CTR])
    print("\nCTR table data is similar")
    assert data['tables']['EXP'][common_columns_EXP].equals(data_R_EXP[common_columns_EXP])
    print("EXP table data is similar")
    assert df1_sorted.equals(df2_sorted)
    print("EXP_x_CTR table data is similar")

def test_check_rankings_data():
    print("\n")
    rankings_list = ['CTR', 'EXP', 'EXP_x_CTR', 'EXP_x_CTR_filtered', 'CTR_ggi', 'EXP_ggi', 'EXP_x_CTR_ggi']
    for ranking in rankings_list:
        # Load table data from R
        ranking_R = pd.read_csv(f'R_data/ranking_{ranking}.csv')
        ranking_R = ranking_R.sort_values(by='nodes').reset_index(drop=True)
        
        ranking_py = data['rankings'][ranking].sort_values(by='nodes').reset_index(drop=True)
        ranking_py['Mediator'] = ranking_py['Mediator'].astype(int)

        assert ranking_py['nodes'].equals(ranking_R['nodes'])
        assert ranking_py['Listener'].equals(ranking_R['Listener'])
        assert ranking_py['Influencer'].equals(ranking_R['Influencer'])
        assert ranking_py['Mediator'].equals(ranking_R['Mediator'])
        assert (abs(ranking_py['Pagerank'] - ranking_R['Pagerank']) < 0.01).all()
        print(f"{ranking} ranking data is similar")

def test_check_graphs_data():
    print("\n")
    graphs_list = ['CTR', 'EXP', 'EXP_x_CTR',  'EXP_x_CTR_filtered', 'CTR_ggi', 'EXP_ggi', 'EXP_x_CTR_ggi']
    for i, graph in enumerate(graphs_list):
        # Load table data from R
        graph_R = pd.read_csv(f'R_data/graph_{graph}.csv')
        graph_R_df = pd.DataFrame({
            "edge": graph_R['from'] + '-' + graph_R['to'],
            "LRScore": graph_R['LRScore'],
        })

        graph_py = data['graphs'][graph] if i < 4  else data['graphs_ggi'][str(graph)[:-4]]
        edges = [edge[0] + '-' + edge[1] for edge in graph_py.edges()]
        lrscores = [edge[2]['LRScore'] for edge in graph_py.edges(data=True)]
        graph_py_df = pd.DataFrame({
            "edge": edges,
            "LRScore": lrscores,
        })

        graph_R_df = graph_R_df.sort_values(by='edge').reset_index(drop=True)
        graph_py_df = graph_py_df.sort_values(by='edge').reset_index(drop=True)

        assert graph_py_df['edge'].equals(graph_R_df['edge'])
        assert (abs(graph_py_df['LRScore'] - graph_R_df['LRScore']) < 0.01).all()

        print(f"{graph} graph data is similar")

def test_check_stats_data():
    stats_R = pd.read_csv('R_data/stat_EXP_x_CTR.csv')
    stats_R_df = pd.DataFrame({
        "cellpair": stats_R['columns_name'],
        "p_value": stats_R['p'],
        "lodds": stats_R['lodds'],
    })
    stats_R_df = stats_R_df.sort_values(by='cellpair').reset_index(drop=True)

    stats_py = data['stats']['EXP_x_CTR']
    stats_py_df = stats_py.sort_values(by='cellpair').reset_index(drop=True)

    assert stats_py_df['cellpair'].equals(stats_R_df['cellpair'])
    assert (abs(stats_py_df['p_value'] - stats_R_df['p_value']) < 0.01).all()
    assert (abs(stats_py_df['lodds'] - stats_R_df['lodds']) < 0.01).all()

    print("\nEXP_x_CTR stats data is similar")

    for key in ['EXP', 'EXP_x_CTR']:
        stats_mannu_R = pd.read_csv(f'R_data/stat_{key}_x_CTR:MannU.csv')
        stats_mannu_R_df = pd.DataFrame({
            "cellpair": stats_mannu_R['cellpair'],
            "p_value": stats_mannu_R['p'],
        })
        stats_mannu_R_df = stats_mannu_R_df.sort_values(by='cellpair').reset_index(drop=True)

        stats_mannu_py = data['stats'][f'{key}_x_CTR:MannU']
        stats_mannu_py_df = stats_mannu_py.sort_values(by='cellpair').reset_index(drop=True)

        assert stats_mannu_py_df['cellpair'].equals(stats_mannu_R_df['cellpair'])
        assert (abs(stats_mannu_py_df['p_value'] - stats_mannu_R_df['p_value']) < 0.01).all()

        print(f"{key}_x_CTR:MannU stats data is similar")