"""
Comprehensive tests to improve coverage for pyCrossTalkeR.

Covers uncovered paths in:
  - pycrosstalker/tools/utils.py
  - pycrosstalker/tools/Single_Condition.py
  - pycrosstalker/tools/Comparative_condition.py
  - pycrosstalker/plots/plot.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from pycrosstalker import plots as ctpl
from pycrosstalker import tools as cttl
from pycrosstalker.tools.utils import (
    add_node_type,
    comparative_med,
    comparative_pagerank,
    create_ordered_circular_layout,
    fisher_test_cci,
    from_liana,
    get_clustered_node_order,
    graph_from_storable_json,
    graph_to_storable_json,
    mannwhitneyu_test_cci,
    ranking_net,
)

# ---------------------------------------------------------------------------
# Shared constants and helpers (also defined in conftest.py for fixtures)
# ---------------------------------------------------------------------------

SEL_COLUMNS = ['source', 'target', 'gene_A', 'gene_B', 'type_gene_A', 'type_gene_B', 'MeanLR']


def make_lr_df(seed=0):
    """Create a minimal synthetic LR interaction DataFrame.

    Uses 3 cell types and 3 LR pairs so PCA(n_components=2) always succeeds
    (3 nodes >= 2 required components).
    """
    np.random.seed(seed)
    cell_types = ['CellA', 'CellB', 'CellC']
    lr_pairs = [('GENE_L1', 'GENE_R1'), ('GENE_L2', 'GENE_R2'), ('GENE_L3', 'GENE_R3')]
    rows = []
    for src in cell_types:
        for tgt in cell_types:
            for gl, gr in lr_pairs:
                rows.append({
                    'source': src,
                    'target': tgt,
                    'gene_A': gl,
                    'gene_B': gr,
                    'type_gene_A': 'Ligand',
                    'type_gene_B': 'Receptor',
                    'MeanLR': abs(float(np.random.normal(1.5, 0.5))) + 0.1,
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Local helpers
# ---------------------------------------------------------------------------

def _make_analysed_pair():
    """Return an AnnData after read_lr_single_condition + create_diff_table."""
    adata = AnnData()
    adata.uns['pycrosstalker'] = {
        'path': {'CTR': make_lr_df(0), 'EXP': make_lr_df(1)}
    }
    adata = cttl.read_lr_single_condition(adata, sel_columns=SEL_COLUMNS)
    adata = cttl.create_diff_table(adata, "/tmp/", comparison=None)
    return adata


def _close():
    plt.close('all')


# ===========================================================================
# utils.py – add_node_type
# ===========================================================================

class TestAddNodeType:
    def test_ligand_receptor(self):
        df = pd.DataFrame({
            'gene_A': ['GX'], 'gene_B': ['GY'],
            'type_gene_A': ['Ligand'], 'type_gene_B': ['Receptor'],
        })
        out = add_node_type(df)
        assert out['gene_A'].iloc[0] == 'GX|L'
        assert out['gene_B'].iloc[0] == 'GY|R'

    def test_receptor_ligand(self):
        df = pd.DataFrame({
            'gene_A': ['GX'], 'gene_B': ['GY'],
            'type_gene_A': ['Receptor'], 'type_gene_B': ['Ligand'],
        })
        out = add_node_type(df)
        assert out['gene_A'].iloc[0] == 'GX|R'
        assert out['gene_B'].iloc[0] == 'GY|L'

    def test_transcription_factor(self):
        df = pd.DataFrame({
            'gene_A': ['TF1'], 'gene_B': ['TF2'],
            'type_gene_A': ['Transcription Factor'],
            'type_gene_B': ['Transcription Factor'],
        })
        out = add_node_type(df)
        assert out['gene_A'].iloc[0] == 'TF1|TF'
        assert out['gene_B'].iloc[0] == 'TF2|TF'


# ===========================================================================
# utils.py – graph_to_storable_json / graph_from_storable_json
# ===========================================================================

class TestGraphSerialization:
    def _sample_graph(self):
        G = nx.DiGraph()
        G.add_edge('CellA', 'CellB', LRScore=1.5, weight=1.5)
        G.add_edge('CellB', 'CellC', LRScore=2.0, weight=2.0)
        return G

    def test_to_json_returns_string(self):
        s = graph_to_storable_json(self._sample_graph())
        assert isinstance(s, str)

    def test_json_is_parseable(self):
        s = graph_to_storable_json(self._sample_graph())
        parsed = json.loads(s)
        assert 'nodes' in parsed
        assert 'edges' in parsed

    def test_roundtrip_nodes(self):
        G = self._sample_graph()
        G2 = graph_from_storable_json(graph_to_storable_json(G))
        assert set(str(n) for n in G2.nodes()) == set(G.nodes())

    def test_roundtrip_edges(self):
        G = self._sample_graph()
        G2 = graph_from_storable_json(graph_to_storable_json(G))
        orig_edges = {(str(u), str(v)) for u, v in G.edges()}
        new_edges = {(str(u), str(v)) for u, v in G2.edges()}
        assert orig_edges == new_edges

    def test_edge_attrs_preserved(self):
        G = nx.DiGraph()
        G.add_edge('X', 'Y', LRScore=3.14, weight=3.14)
        G2 = graph_from_storable_json(graph_to_storable_json(G))
        assert pytest.approx(G2['X']['Y']['LRScore'], abs=1e-5) == 3.14


# ===========================================================================
# utils.py – get_clustered_node_order / create_ordered_circular_layout
# ===========================================================================

class TestLayoutUtils:
    def _edge_df(self):
        return pd.DataFrame({
            'source': ['A', 'B', 'C', 'A'],
            'target': ['B', 'C', 'A', 'C'],
            'LRScore': [1.0, 2.0, 0.5, 1.5],
        })

    def test_get_clustered_node_order_returns_all_nodes(self):
        order = get_clustered_node_order(self._edge_df())
        assert isinstance(order, list)
        assert set(order) == {'A', 'B', 'C'}
        assert len(order) == 3

    def test_get_clustered_node_order_no_duplicates(self):
        order = get_clustered_node_order(self._edge_df())
        assert len(order) == len(set(order))

    def test_create_ordered_circular_layout_count(self):
        layout = create_ordered_circular_layout(['A', 'B', 'C', 'D'])
        assert len(layout) == 4

    def test_create_ordered_circular_layout_on_unit_circle(self):
        layout = create_ordered_circular_layout(['A', 'B', 'C'])
        for node, (x, y) in layout.items():
            assert abs(x**2 + y**2 - 1.0) < 1e-9, f"Node {node} not on unit circle"


# ===========================================================================
# utils.py – ranking_net (unsigned and signed modes)
# ===========================================================================

class TestRankingNet:
    def _simple_digraph(self):
        G = nx.DiGraph()
        G.add_edge('A', 'B', weight=1.0)
        G.add_edge('B', 'C', weight=2.0)
        G.add_edge('A', 'C', weight=0.5)
        return G

    def test_unsigned_columns(self):
        result = ranking_net(self._simple_digraph(), mode=True)
        assert set(result.columns) >= {'nodes', 'Listener', 'Influencer', 'Mediator'}

    def test_unsigned_all_nodes_present(self):
        result = ranking_net(self._simple_digraph(), mode=True)
        assert set(result['nodes']) == {'A', 'B', 'C'}

    def test_unsigned_no_nan(self):
        result = ranking_net(self._simple_digraph(), mode=True)
        assert not result[['Listener', 'Influencer', 'Mediator']].isna().any().any()

    def test_signed_columns(self):
        G = nx.DiGraph()
        G.add_edge('A', 'B', weight=1.0)
        G.add_edge('B', 'C', weight=-2.0)
        result = ranking_net(G, mode=False)
        assert set(result.columns) >= {'nodes', 'Listener', 'Influencer'}

    def test_signed_no_mediator_column(self):
        G = nx.DiGraph()
        G.add_edge('A', 'B', weight=1.0)
        result = ranking_net(G, mode=False)
        assert 'Mediator' not in result.columns

    def test_signed_all_nodes_present(self):
        G = nx.DiGraph()
        G.add_edge('A', 'B', weight=1.0)
        G.add_edge('B', 'C', weight=-2.0)
        result = ranking_net(G, mode=False)
        assert set(result['nodes']) == {'A', 'B', 'C'}


# ===========================================================================
# utils.py – comparative_pagerank / comparative_med
# ===========================================================================

class TestComparativeRankingHelpers:
    def _rankings_fixture(self):
        """Minimal rankings dict with CTR and EXP entries."""
        ctr = pd.DataFrame({'nodes': ['X', 'Y'], 'Pagerank': [0.6, 0.4], 'Mediator': [5.0, 3.0]})
        exp = pd.DataFrame({'nodes': ['X', 'Y'], 'Pagerank': [0.7, 0.3], 'Mediator': [6.0, 2.0]})
        return {'CTR': ctr, 'EXP': exp}

    def _curr_rkg(self):
        return pd.DataFrame({'nodes': ['X', 'Y'], 'Listener': [1, 2], 'Influencer': [2, 1]})

    def test_comparative_pagerank_adds_column(self):
        rankings = self._rankings_fixture()
        curr = self._curr_rkg()
        result = comparative_pagerank(rankings, "graphs", "EXP_x_CTR", curr)
        assert 'Pagerank' in result.columns

    def test_comparative_pagerank_filtered_name(self):
        rankings = self._rankings_fixture()
        curr = self._curr_rkg()
        result = comparative_pagerank(rankings, "graphs", "EXP_x_CTR_filtered", curr)
        assert 'Pagerank' in result.columns

    def test_comparative_med_adds_column(self):
        rankings = self._rankings_fixture()
        curr = self._curr_rkg()
        result = comparative_med(rankings, "graphs", "EXP_x_CTR", curr)
        assert 'Mediator' in result.columns

    def test_comparative_med_filtered_name(self):
        rankings = self._rankings_fixture()
        curr = self._curr_rkg()
        result = comparative_med(rankings, "graphs", "EXP_x_CTR_filtered", curr)
        assert 'Mediator' in result.columns


# ===========================================================================
# utils.py – fisher_test_cci with explicit comparison parameter
# ===========================================================================

class TestFisherTestExplicitComparison:
    def test_explicit_comparison_adds_stats(self):
        adata = _make_analysed_pair()
        adata = fisher_test_cci(adata, 'LRScore', '/tmp/', comparison=[('EXP', 'CTR')])
        stats = adata.uns['pycrosstalker']['results']['stats']
        assert 'EXP_x_CTR' in stats

    def test_explicit_comparison_result_has_expected_columns(self):
        adata = _make_analysed_pair()
        adata = fisher_test_cci(adata, 'LRScore', '/tmp/', comparison=[('EXP', 'CTR')])
        result = adata.uns['pycrosstalker']['results']['stats']['EXP_x_CTR']
        assert {'cellpair', 'p_value', 'lodds'}.issubset(result.columns)

    def test_explicit_comparison_p_values_in_range(self):
        adata = _make_analysed_pair()
        adata = fisher_test_cci(adata, 'LRScore', '/tmp/', comparison=[('EXP', 'CTR')])
        p_vals = adata.uns['pycrosstalker']['results']['stats']['EXP_x_CTR']['p_value']
        assert ((p_vals >= 0) & (p_vals <= 1)).all()


# ===========================================================================
# utils.py – mannwhitneyu_test_cci with explicit comparison parameter
# ===========================================================================

class TestMannWhitneyExplicitComparison:
    def test_explicit_comparison_adds_stats(self):
        adata = _make_analysed_pair()
        adata = mannwhitneyu_test_cci(adata, 'LRScore', '/tmp/', comparison=[('EXP', 'CTR')])
        stats = adata.uns['pycrosstalker']['results']['stats']
        assert 'EXP_x_CTR:MannU' in stats

    def test_explicit_comparison_result_columns(self):
        adata = _make_analysed_pair()
        adata = mannwhitneyu_test_cci(adata, 'LRScore', '/tmp/', comparison=[('EXP', 'CTR')])
        result = adata.uns['pycrosstalker']['results']['stats']['EXP_x_CTR:MannU']
        assert {'cellpair', 'statistic', 'p_value', 'lfc'}.issubset(result.columns)


# ===========================================================================
# utils.py – from_liana
# ===========================================================================

class TestFromLiana:
    def _base_df(self):
        return pd.DataFrame({
            'ligand': ['GENE_L1', 'GENE_L2'],
            'receptor_complex': ['GENE_R1', 'GENE_R2'],
            'source': ['CellA', 'CellB'],
            'target': ['CellB', 'CellA'],
            'cellphone_pvals': [0.01, 0.03],
            'lr_means': [1.5, 2.0],
        })

    def test_dataframe_with_label_creates_path(self):
        df = self._base_df()
        df['label'] = 'condition1_lr_data'
        adata = AnnData()
        adata.uns['liana'] = df
        result = from_liana(adata)
        assert 'pycrosstalker' in result.uns
        assert len(result.uns['pycrosstalker']['path']) > 0

    def test_dataframe_pval_filter_applied(self):
        df = self._base_df()
        df['label'] = 'condition1_lr_data'
        # Both p_values <= 0.05, so both should pass
        adata = AnnData()
        adata.uns['liana'] = df
        result = from_liana(adata, pval_filter=True)
        key = list(result.uns['pycrosstalker']['path'].keys())[0]
        assert len(result.uns['pycrosstalker']['path'][key]) == 2

    def test_dataframe_no_pval_filter(self):
        df = self._base_df()
        df['label'] = 'condition1_lr_data'
        df['cellphone_pvals'] = [0.5, 0.8]  # both would be filtered out
        adata = AnnData()
        adata.uns['liana'] = df
        result = from_liana(adata, pval_filter=False)
        key = list(result.uns['pycrosstalker']['path'].keys())[0]
        assert len(result.uns['pycrosstalker']['path'][key]) == 2

    def test_dict_input_creates_path(self):
        adata = AnnData()
        adata.uns['liana'] = {'condition1_lr_data': self._base_df()}
        result = from_liana(adata, pval_filter=False)
        assert 'pycrosstalker' in result.uns
        assert 'condition1_lr_data' in result.uns['pycrosstalker']['path']

    def test_dict_correct_columns(self):
        adata = AnnData()
        adata.uns['liana'] = {'condition1_lr_data': self._base_df()}
        result = from_liana(adata, pval_filter=False)
        df = result.uns['pycrosstalker']['path']['condition1_lr_data']
        expected = {'source', 'target', 'type_gene_A', 'type_gene_B', 'gene_A', 'gene_B', 'MeanLR'}
        assert expected.issubset(df.columns)

    def test_compute_means(self):
        df = self._base_df()
        df['label'] = 'condition1_lr_data'
        df['ligand_means'] = [0.8, 1.2]
        df['receptor_means'] = [1.1, 0.9]
        adata = AnnData()
        adata.uns['liana'] = df
        result = from_liana(adata, compute_means=True, pval_filter=False)
        assert len(result.uns['pycrosstalker']['path']) > 0
        key = list(result.uns['pycrosstalker']['path'].keys())[0]
        assert 'MeanLR' in result.uns['pycrosstalker']['path'][key].columns


# ===========================================================================
# Single_Condition.py – dict input, invalid input, receptor-first ordering
# ===========================================================================

class TestReadLRSingleCondition:
    def test_dict_input_from_csv(self, tmp_path):
        csv_path = tmp_path / "CTR_LR.csv"
        make_lr_df(0).to_csv(csv_path, index=False)
        result = cttl.read_lr_single_condition({'CTR': str(csv_path)}, sel_columns=SEL_COLUMNS)
        assert isinstance(result, AnnData)
        assert 'CTR' in result.uns['pycrosstalker']['results']['graphs']

    def test_dict_input_creates_correct_structure(self, tmp_path):
        csv_path = tmp_path / "CTR_LR.csv"
        make_lr_df(0).to_csv(csv_path, index=False)
        result = cttl.read_lr_single_condition({'CTR': str(csv_path)}, sel_columns=SEL_COLUMNS)
        keys = {'graphs', 'graphs_ggi', 'tables', 'colors', 'coords', 'rankings', 'pca', 'stats'}
        assert keys.issubset(result.uns['pycrosstalker']['results'].keys())

    def test_invalid_input_type_raises(self):
        with pytest.raises(ValueError, match="Input parameter"):
            cttl.read_lr_single_condition(42, sel_columns=SEL_COLUMNS)

    def test_invalid_dict_value_raises(self):
        with pytest.raises(ValueError, match="Issue with input paths"):
            cttl.read_lr_single_condition({'CTR': 42}, sel_columns=SEL_COLUMNS)

    def test_receptor_first_ordering(self):
        """Covers the else-branch when type_gene_A == 'Receptor' (not 'Ligand')."""
        df = pd.DataFrame({
            'source': ['CellA', 'CellA', 'CellB'],
            'target': ['CellB', 'CellB', 'CellA'],
            'gene_A': ['GENE_R1', 'GENE_R2', 'GENE_R1'],
            'gene_B': ['GENE_L1', 'GENE_L2', 'GENE_L1'],
            'type_gene_A': ['Receptor', 'Receptor', 'Receptor'],
            'type_gene_B': ['Ligand', 'Ligand', 'Ligand'],
            'MeanLR': [1.5, 2.0, 0.8],
        })
        adata = AnnData()
        adata.uns['pycrosstalker'] = {'path': {'CTR': df}}
        result = cttl.read_lr_single_condition(adata, sel_columns=SEL_COLUMNS)
        assert 'CTR' in result.uns['pycrosstalker']['results']['graphs']
        # ligpair/recpair should be correctly formed even with swapped gene order
        assert 'ligpair' in result.uns['pycrosstalker']['results']['tables']['CTR'].columns


# ===========================================================================
# Comparative_condition.py – explicit comparison parameter
# ===========================================================================

class TestCreateDiffTableExplicit:
    def test_explicit_comparison_creates_diff_graph(self):
        adata = AnnData()
        adata.uns['pycrosstalker'] = {'path': {'CTR': make_lr_df(0), 'EXP': make_lr_df(1)}}
        adata = cttl.read_lr_single_condition(adata, sel_columns=SEL_COLUMNS)
        adata = cttl.create_diff_table(adata, "/tmp/", comparison=[('EXP', 'CTR')])
        results = adata.uns['pycrosstalker']['results']
        assert 'EXP_x_CTR' in results['graphs']
        assert 'EXP_x_CTR' in results['tables']
        assert 'EXP_x_CTR' in results['graphs_ggi']

    def test_explicit_comparison_lrscore_column_present(self):
        adata = AnnData()
        adata.uns['pycrosstalker'] = {'path': {'CTR': make_lr_df(0), 'EXP': make_lr_df(1)}}
        adata = cttl.read_lr_single_condition(adata, sel_columns=SEL_COLUMNS)
        adata = cttl.create_diff_table(adata, "/tmp/", comparison=[('EXP', 'CTR')])
        graph_df = adata.uns['pycrosstalker']['results']['graphs']['EXP_x_CTR']
        assert 'LRScore' in graph_df.columns

    def test_explicit_multiple_comparisons(self):
        adata = AnnData()
        adata.uns['pycrosstalker'] = {
            'path': {'CTR': make_lr_df(0), 'EXP1': make_lr_df(1), 'EXP2': make_lr_df(2)}
        }
        adata = cttl.read_lr_single_condition(adata, sel_columns=SEL_COLUMNS)
        adata = cttl.create_diff_table(
            adata, "/tmp/", comparison=[('EXP1', 'CTR'), ('EXP2', 'CTR')]
        )
        results = adata.uns['pycrosstalker']['results']
        assert 'EXP1_x_CTR' in results['graphs']
        assert 'EXP2_x_CTR' in results['graphs']


# ===========================================================================
# plots/plot.py – plot_cci additional branches
# ===========================================================================

class TestPlotCCIBranches:
    def test_plot_cci_log_true(self, analysed_adata):
        """Covers the log=True branch (line ~104)."""
        data = analysed_adata.uns['pycrosstalker']['results']
        fig, ax = ctpl.plot.plot_cci(
            graph=data['graphs']['CTR'],
            colors=data['colors'],
            plt_name='Log Mode',
            coords=data['coords'],
            pg=list(data['rankings']['CTR']['Pagerank']),
            log=True,
            return_figure=True,
        )
        assert fig is not None
        _close()

    def test_plot_cci_single_coord(self):
        """Covers the single-node coordinate branch (line ~85)."""
        graph_df = pd.DataFrame({
            'source': ['CellA'],
            'target': ['CellA'],
            'LRScore': [1.0],
            'freq': [0.5],
            'weight': [1.0],
            'inter': [0.5],
        })
        colors = {'CellA': '#ff0000'}
        coords = {'CellA': (0.0, 0.0)}
        fig, ax = ctpl.plot.plot_cci(
            graph=graph_df,
            colors=colors,
            plt_name='Single Node',
            coords=coords,
            pg=[0.5],
            return_figure=True,
        )
        assert fig is not None
        _close()

    def test_plot_cci_with_emax(self, analysed_adata):
        """Covers emax-specified branch."""
        data = analysed_adata.uns['pycrosstalker']['results']
        fig, ax = ctpl.plot.plot_cci(
            graph=data['graphs']['CTR'],
            colors=data['colors'],
            plt_name='emax set',
            coords=data['coords'],
            pg=list(data['rankings']['CTR']['Pagerank']),
            emax=5.0,
            return_figure=True,
        )
        assert fig is not None
        _close()


# ===========================================================================
# plots/plot.py – plot_bar_rankings
# ===========================================================================

class TestPlotBarRankings:
    def test_comparative_ranking_runs(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank'
        )
        _close()

    def test_comparative_ranking_mode_cgi(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank', mode='cgi'
        )
        _close()

    def test_non_comparative_returns_none(self, analysed_adata):
        """plot_bar_rankings with no _x_ table silently returns None."""
        result = ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='CTR', ranking='Pagerank'
        )
        assert result is None
        _close()

    def test_filter_sign_neg(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank', filter_sign='neg'
        )
        _close()


# ===========================================================================
# plots/plot.py – plot_clustermap / plot_graph_clustermap
# ===========================================================================

class TestPlotClustermaps:
    def test_plot_clustermap_runs(self, analysed_adata):
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR']
        ctpl.plot.plot_clustermap(table, title='CTR Clustermap')
        _close()

    def test_plot_graph_clustermap_runs(self, analysed_adata):
        graph = analysed_adata.uns['pycrosstalker']['results']['graphs']['CTR']
        ctpl.plot.plot_graph_clustermap(graph, title='CTR Graph Clustermap')
        _close()

    def test_plot_graph_clustermap_differential(self, analysed_adata):
        """Covers the comparative (signed) graph case with negative values."""
        graph = analysed_adata.uns['pycrosstalker']['results']['graphs']['EXP_x_CTR']
        ctpl.plot.plot_graph_clustermap(graph, title='EXP_x_CTR Graph Clustermap')
        _close()


# ===========================================================================
# plots/plot.py – plot_volcane
# ===========================================================================

class TestPlotVolcane:
    def test_fisher_volcane_runs(self, analysed_adata):
        stats = analysed_adata.uns['pycrosstalker']['results']['stats']['EXP_x_CTR'].copy()
        ctpl.plot.plot_volcane(stats, method='fisher')
        _close()

    def test_mannwhitney_volcane_runs(self, analysed_adata):
        stats = analysed_adata.uns['pycrosstalker']['results']['stats'].get('EXP_x_CTR:MannU')
        if stats is not None and len(stats) > 0:
            ctpl.plot.plot_volcane(stats.copy(), method='mannwhitneyu')
        _close()


# ===========================================================================
# plots/plot.py – plot_pca_LR_comparative
# ===========================================================================

class TestPlotPCAComparative:
    def test_pca_non_ggi_returns_plot(self, analysed_adata):
        data = analysed_adata.uns['pycrosstalker']['results']
        result = ctpl.plot.plot_pca_LR_comparative(data, 'CTR', ret=True, ggi=False)
        assert result is not None
        _close()

    def test_pca_ggi_all_gene_types(self, analysed_adata):
        data = analysed_adata.uns['pycrosstalker']['results']
        ctpl.plot.plot_pca_LR_comparative(data, 'CTR_ggi', ret=False, ggi=True, gene_types='all')
        _close()

    def test_pca_ggi_lr_gene_types(self, analysed_adata):
        data = analysed_adata.uns['pycrosstalker']['results']
        ctpl.plot.plot_pca_LR_comparative(data, 'CTR_ggi', ret=False, ggi=True, gene_types='LR')
        _close()

    def test_pca_ggi_include_tf(self, analysed_adata):
        import copy
        # Deep-copy the rankings/pca entries so in-place PC1/PC2 writes from earlier
        # tests do not corrupt the shared session fixture for this or later tests.
        data = analysed_adata.uns['pycrosstalker']['results']
        data_copy = dict(data)
        data_copy['rankings'] = {k: v.copy() for k, v in data['rankings'].items()}
        ctpl.plot.plot_pca_LR_comparative(
            data_copy, 'CTR_ggi', ret=False, ggi=True, include_tf=True, gene_types='all'
        )
        _close()

    def test_pca_ggi_returns_plot_when_ret_true(self, analysed_adata):
        data = analysed_adata.uns['pycrosstalker']['results']
        result = ctpl.plot.plot_pca_LR_comparative(data, 'CTR_ggi', ret=True, ggi=True)
        assert result is not None
        _close()


# ===========================================================================
# plots/plot.py – plot_sankey / gen_sankey
# ===========================================================================

class TestPlotSankey:
    def test_no_target_uses_all_data(self, analysed_adata):
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR'].copy()
        ctpl.plot.plot_sankey(table, target=None, plt_name='All interactions')
        _close()

    def test_with_ligand_target(self, analysed_adata):
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR'].copy()
        ligand_genes = table[table['type_gene_A'] == 'Ligand']['gene_A'].dropna().unique()
        if len(ligand_genes) > 0:
            ctpl.plot.plot_sankey(table.copy(), target=ligand_genes[0], plt_name='Ligand target')
        _close()

    def test_with_receptor_target(self, analysed_adata):
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR'].copy()
        receptor_genes = table[table['type_gene_B'] == 'Receptor']['gene_B'].dropna().unique()
        if len(receptor_genes) > 0:
            ctpl.plot.plot_sankey(table.copy(), target=receptor_genes[0], plt_name='Receptor target')
        _close()

    def test_with_ligand_cluster_filter(self, analysed_adata):
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR'].copy()
        ctpl.plot.plot_sankey(
            table, target=None, ligand_cluster=['CellA'], plt_name='Ligand cluster'
        )
        _close()

    def test_with_receptor_cluster_filter(self, analysed_adata):
        """Covers the receptor_cluster is not None branch (line ~433)."""
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR'].copy()
        ctpl.plot.plot_sankey(
            table, target=None, receptor_cluster=['CellB'], plt_name='Receptor cluster'
        )
        _close()

    def test_with_missing_target_prints_not_found(self, analysed_adata, capsys):
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR'].copy()
        ctpl.plot.plot_sankey(table, target='NONEXISTENT_GENE_XYZ', plt_name='Missing')
        captured = capsys.readouterr()
        assert 'Not Found' in captured.out
        _close()


# ===========================================================================
# utils.py – remaining gaps
# ===========================================================================

class TestGraphSerializationNumpyScalars:
    def test_numpy_scalar_edge_attrs_serialized(self):
        """Covers the isinstance(v, np.generic) branch in graph_to_storable_json (line ~541)."""
        G = nx.DiGraph()
        G.add_edge('A', 'B', LRScore=np.float64(2.5), weight=np.float64(2.5))
        json_str = graph_to_storable_json(G)
        # Must round-trip without TypeError (np.float64 is not JSON-serializable natively)
        parsed = json.loads(json_str)
        assert 'edges' in parsed


class TestFromLianaRemainingBranches:
    def _base_df(self):
        return pd.DataFrame({
            'ligand': ['GENE_L1', 'GENE_L2'],
            'receptor_complex': ['GENE_R1', 'GENE_R2'],
            'source': ['CellA', 'CellB'],
            'target': ['CellB', 'CellA'],
            'cellphone_pvals': [0.01, 0.03],
            'lr_means': [1.5, 2.0],
            'ligand_means': [0.8, 1.2],
            'receptor_means': [1.1, 0.9],
        })

    def test_dict_compute_means(self):
        """Covers compute_means=True in the dict path (line ~737)."""
        adata = AnnData()
        adata.uns['liana'] = {'condition1_lr_data': self._base_df()}
        result = from_liana(adata, compute_means=True, pval_filter=False)
        key = list(result.uns['pycrosstalker']['path'].keys())[0]
        assert 'MeanLR' in result.uns['pycrosstalker']['path'][key].columns

    def test_dict_pval_filter_true(self):
        """Covers pval_filter=True in the dict path (line ~740)."""
        df = self._base_df().copy()
        df['cellphone_pvals'] = [0.01, 0.8]  # only first row passes
        adata = AnnData()
        adata.uns['liana'] = {'condition1': df}
        result = from_liana(adata, pval_filter=True)
        key = list(result.uns['pycrosstalker']['path'].keys())[0]
        assert len(result.uns['pycrosstalker']['path'][key]) == 1


# ===========================================================================
# plots/plot.py – plot_bar_rankings remaining branches
# ===========================================================================

class TestPlotBarRankingsTypeFiler:
    """Cover the 'type' parameter branches in plot_bar_rankings (lines 314-326)."""

    def test_type_single_char_L(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank', type='L'
        )
        _close()

    def test_type_two_char_lr(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank', type='LR'
        )
        _close()

    def test_type_two_char_tf(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank', type='TF'
        )
        _close()

    def test_type_three_char_rtf(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank', type='RTF'
        )
        _close()

    def test_type_three_char_ltf(self, analysed_adata):
        ctpl.plot.plot_bar_rankings(
            analysed_adata, table_name='EXP_x_CTR', ranking='Pagerank', type='LTF'
        )
        _close()

    def test_filter_sign_pos_empty_returns_message(self, analysed_adata):
        """Covers filter_sign='pos' (line 337) and empty check (line 343)."""
        # Use a ranking column where all values for the filtered set are negative
        # by filtering for a type that doesn't exist in the data, guaranteeing empty
        result = ctpl.plot.plot_bar_rankings(
            analysed_adata,
            table_name='EXP_x_CTR',
            ranking='Pagerank',
            type='TF',       # No TF nodes → empty table
            filter_sign='pos',
        )
        assert result == "No entries with provided Filters."
        _close()


# ===========================================================================
# plots/plot.py – plot_volcane with annotated (red) points
# ===========================================================================

class TestPlotVolcaneAnnotated:
    def test_volcane_with_red_points(self):
        """Covers the text annotation branch when points are significant (line ~760)."""
        # Construct a stats DataFrame where at least one cell pair is 'red'
        # (|lfc| > fc_threshold AND p_value < p_threshold)
        df = pd.DataFrame({
            'cellpair': ['CellA@CellB', 'CellA@CellA', 'CellB@CellB'],
            'p_value':  [0.001, 0.5, 0.8],
            'lodds':    [3.0, 0.1, -0.2],  # first entry has |lodds|>1 and p<0.05
        })
        ctpl.plot.plot_volcane(df.copy(), method='fisher', p_threshold=0.05, fc_threshold=1)
        _close()


# ===========================================================================
# plots/plot.py – gen_sankey2 (alternative Plotly-based Sankey)
# ===========================================================================

class TestGenSankey2:
    def test_gen_sankey2_runs(self, analysed_adata):
        """Directly exercises gen_sankey2 (lines 474-557)."""
        table = analysed_adata.uns['pycrosstalker']['results']['tables']['CTR'].copy()
        # gen_sankey2 requires positive count values; filter to ligand/receptor rows
        df = table[
            (table['type_gene_A'] == 'Ligand') & (table['type_gene_B'] == 'Receptor')
        ].head(5).copy()
        df = df[['source', 'gene_A', 'gene_B', 'target', 'LRScore']].copy()
        df['LRScore'] = df['LRScore'].abs() + 0.1  # ensure positive, non-zero
        ctpl.plot.gen_sankey2(
            df, cat_cols=['source', 'gene_A', 'gene_B', 'target'], value_cols='LRScore',
            title='Test gen_sankey2'
        )
        _close()
