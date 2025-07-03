from pycrosstalker import tools as cttl
from pycrosstalker import plots as ctpl
import pickle

# Load data from pickle file
with open("../tutorials/output/HMI_example/LR_data.pkl", "rb") as f:
    data = pickle.load(f)

def test_data_type():
    print("\nTesting data type is dict")
    assert isinstance(data, dict)

def test_data_keys():
    print("Testing data keys")
    keys = {'graphs', 'graphs_ggi', 'tables', 'colors', 'coords', 'rankings', 'pca', 'stats'}
    assert keys.issubset(data.keys())