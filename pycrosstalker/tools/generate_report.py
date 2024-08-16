import os
import pandas as pd
from tqdm import tqdm
from .Comparative_condition import *
from .Single_Condition import *
from .utils import *


def analise_LR(lrpaths, 
                genes=None, 
                tf_genes=None, 
                out_path=None, 
                sep=',', 
                threshold=0, 
                colors=None, 
                out_file=None, 
                output_fmt="html_document", 
                sel_columns=['source','target','gene_A','gene_B','type_gene_A','type_gene_B','MeanLR'], 
                org='hsa', comparison=None, filtered_net=False):
    
    data = read_lr_single_condition(lrpaths, 
                                    sel_columns, 
                                    out_path, 
                                    sep, 
                                    colors)

    print("Create a Differential Table")
    if len(lrpaths) > 1:
        data = create_diff_table(data, out_path, comparison)

    print("Calculating CCI Ranking")
    data = ranking(data, out_path, sel_columns=sel_columns, slot="graphs")
    print("Calculating GCI Ranking")
    data = ranking(data, out_path, sel_columns=sel_columns, slot="graphs_ggi")
    print("Network Analysis Done")

    with open(os.path.join(out_path, "LR_data.pkl"), "wb") as f:
        pickle.dump(data, f)

    return(data)

# paths = {
#     'CTR': "Test_Data/CTR_LR.csv",
#     'EXP': "Test_Data/EXP_LR.csv"
# }
# output = "output/"
# data = analise_LR(paths, out_path=output, org="hsa")

