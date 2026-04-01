"""
Shared pytest fixtures for pyCrossTalkeR tests.
"""
import pytest
from anndata import AnnData
from pycrosstalker import tools as cttl
from test.helpers import SEL_COLUMNS, make_lr_df  # noqa: F401 – re-exported for tests

__all__ = ["SEL_COLUMNS", "make_lr_df"]


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
