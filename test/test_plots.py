import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
from pycrosstalker import tools as cttl
from pycrosstalker import plots as ctpl

# Load and analyse data once for all tests
with open("rawdata/humanBM.h5ad", "rb") as f:
    paths = {
        'CTR': "rawdata/CTR_LR.csv",
        'EXP': "rawdata/EXP_LR.csv"
    }
    adata = sc.read_h5ad(f)
    adata.uns['pycrosstalker'] = {}
    adata.uns['pycrosstalker']['path'] = {}
    for k, v in paths.items():
        adata.uns['pycrosstalker']['path'][k] = pd.read_csv(v)
    adata = cttl.analise_LR(adata, save=False)
    data = adata.uns['pycrosstalker']['results']

lrobj_tblPCA = {
    'pca': data['pca'],
    'rankings': data['rankings'],
    'tables': data['tables'],
}


def test_plot_cci():
    fig, ax = ctpl.plot.plot_cci(
        graph=data["graphs"]["CTR"],
        colors=data["colors"],
        plt_name='Control',
        coords=data["coords"],
        emax=None,
        leg=False,
        low=0,
        high=0,
        ignore_alpha=False,
        log=False,
        efactor=6,
        vfactor=12,
        pg=data["rankings"]["CTR"]["Pagerank"],
        scale_factor=2.0,
        node_size=2.5,
        return_figure=True
    )
    assert fig is not None
    assert ax is not None
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_cci_exp():
    fig, ax = ctpl.plot.plot_cci(
        graph=data["graphs"]["EXP"],
        colors=data["colors"],
        plt_name='Experimental',
        coords=data["coords"],
        emax=None,
        leg=False,
        low=0,
        high=0,
        ignore_alpha=False,
        log=True,
        efactor=6,
        vfactor=12,
        pg=data["rankings"]["EXP"]["Pagerank"],
        scale_factor=2.0,
        node_size=2.5,
        return_figure=True
    )
    assert fig is not None
    assert ax is not None
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_pca_LR_comparative_ggi():
    result = ctpl.plot.plot_pca_LR_comparative(
        lrobj_tblPCA=lrobj_tblPCA,
        pca_table='CTR_ggi',
        dims=(1, 2),
        ret=True,
        ggi=True,
        include_tf=False,
        gene_types="all"
    )
    assert result is not None
    assert 'CTR_ggi' in result
    plt.close("all")


def test_plot_pca_LR_comparative_no_ggi():
    result = ctpl.plot.plot_pca_LR_comparative(
        lrobj_tblPCA=lrobj_tblPCA,
        pca_table='CTR',
        dims=(1, 2),
        ret=True,
        ggi=False
    )
    assert result is not None
    assert 'CTR' in result
    plt.close("all")


def test_plot_bar_rankings_cci():
    ctpl.plot.plot_bar_rankings(
        annData=adata,
        table_name='EXP_x_CTR',
        ranking='Pagerank',
        type=None,
        filter_sign=None,
        mode='cci',
        top_num=10
    )
    plt.close("all")


def test_plot_bar_rankings_cgi():
    ctpl.plot.plot_bar_rankings(
        annData=adata,
        table_name='EXP_x_CTR_ggi',
        ranking='Pagerank',
        type=None,
        filter_sign=None,
        mode='cgi',
        top_num=10
    )
    plt.close("all")


def test_plot_bar_rankings_with_type_filter():
    result = ctpl.plot.plot_bar_rankings(
        annData=adata,
        table_name='EXP_x_CTR_ggi',
        ranking='Pagerank',
        type='LR',
        filter_sign=None,
        mode='cgi',
        top_num=10
    )
    # May return "No entries" string if no matches, or None if plotted
    assert result is None or isinstance(result, str)
    plt.close("all")


def test_plot_sankey():
    ctpl.plot.plot_sankey(
        lrobj_tbl=data['tables']['CTR'],
        target=None,
        ligand_cluster=None,
        receptor_cluster=None,
        plt_name='CTR Sankey',
        threshold=50
    )
    plt.close("all")


def test_plot_sankey_with_target():
    # Use the first gene_B value from the CTR table as a receptor target
    sample_gene = data['tables']['CTR']['gene_B'].iloc[0]
    ctpl.plot.plot_sankey(
        lrobj_tbl=data['tables']['CTR'],
        target=sample_gene,
        ligand_cluster=None,
        receptor_cluster=None,
        plt_name='Sankey with target',
        threshold=50
    )
    plt.close("all")


def test_plot_volcane_fisher():
    stats_df = data['stats']['EXP_x_CTR'].copy()
    ctpl.plot.plot_volcane(
        df=stats_df,
        method='fisher',
        p_threshold=0.05,
        fc_threshold=1,
        figsize=(8, 6),
        annot=True,
        title='Volcano Plot (Fisher)'
    )
    plt.close("all")


def test_plot_clustermap():
    ctpl.plot.plot_clustermap(
        data=data['tables']['CTR'],
        title='CTR Clustermap',
        annot=True,
        return_figure=False
    )
    plt.close("all")


def test_plot_graph_clustermap():
    ctpl.plot.plot_graph_clustermap(
        graph=data['graphs']['CTR'],
        weight='LRScore',
        title='CTR Graph Clustermap',
        annot=True
    )
    plt.close("all")


def test_plot_graph_clustermap_exp_x_ctr():
    ctpl.plot.plot_graph_clustermap(
        graph=data['graphs']['EXP_x_CTR'],
        weight='LRScore',
        title='EXP_x_CTR Graph Clustermap',
        annot=True
    )
    plt.close("all")
