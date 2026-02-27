import pandas as pd
import scanpy as sc
from pycrosstalker import tools as cttl
from pycrosstalker import plots as ctpl
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")  # headless backend

# Load data from AnnData file
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


def test_plot_cci():
    fig, ax = ctpl.plot.plot_cci(graph=adata.uns['pycrosstalker']['results']["graphs"]["CTR"],
        colors=adata.uns['pycrosstalker']['results']["colors"],
        plt_name='Control',
        coords=adata.uns['pycrosstalker']['results']["coords"],
        emax=None,
        leg=False,
        low=0,
        high=0,
        ignore_alpha=False,
        log=False,
        efactor=6,
        vfactor=12,
        pg=adata.uns['pycrosstalker']['results']["rankings"]["CTR"]["Pagerank"],
        scale_factor=2.0,
        node_size= 2.5,
        return_figure=True)
    assert fig is not None
    assert ax is not None
    assert isinstance(ax, plt.Axes)
