import scanpy as sc
from anndata import AnnData 
import importlib

def test_import_package():
    importlib.import_module("pycrosstalker")

# Load data from AnnData file
with open("tutorials/output/Myelofibrosis_example/Myelofibrosis_example_analysed.h5ad", "rb") as f:
    adata = sc.read_h5ad(f)

def test_data_type():
    print("\nTesting data type is AnnData")
    assert isinstance(adata, AnnData)

def test_data_keys():
    print("Testing AnnData attributes exit")
    keys = {'graphs', 'graphs_ggi', 'tables', 'colors', 'coords', 'rankings', 'pca', 'stats'}
    assert keys.issubset(adata.uns['pycrosstalker']['results'].keys())


