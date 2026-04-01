"""
Shared pytest fixtures for pyCrossTalkeR tests.
"""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from pycrosstalker import tools as cttl

SEL_COLUMNS = ['source', 'target', 'gene_A', 'gene_B', 'type_gene_A', 'type_gene_B', 'MeanLR']


def make_lr_df(seed=0):
    """Create a minimal synthetic LR interaction DataFrame.

    Uses 3 cell types and 3 LR pairs so that PCA(n_components=2) always has
    enough samples (3 nodes) and features (Listener, Influencer, Mediator).
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


@pytest.fixture(scope="session")
def analysed_adata():
    """Fully analysed AnnData built from synthetic data (runs once per session)."""
    adata = AnnData()
    adata.uns['pycrosstalker'] = {
        'path': {
            'CTR': make_lr_df(0),
            'EXP': make_lr_df(1),
        }
    }
    return cttl.analise_LR(adata, save=False)
