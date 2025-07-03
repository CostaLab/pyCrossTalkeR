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
    
    """
    Core engine to generate report. Here we perform all the computation related to pyCrossTalkeR

    Parameters
    ----------
    lrpaths :
        Paths of single condition LR data
    genes :
        list of genes to be considered in the sankey plots
    out_path :
        output directory path
    sep :
        character used on csv
    threshold :
        percentage of edges to be pruned
    colors :
        celltypes colorscheme
    out_file :
        output file names
    output_fmt :
        rmarkdown render output format parameter
    sel_columns :
        columns from data
    
    Returns
    -------
    Rmarkdown report all objects from each step
    
    """
    
    data = read_lr_single_condition(lrpaths, 
                                    sel_columns, 
                                    out_path, 
                                    sep, 
                                    colors)

    print("Create a Differential Table")
    if len(lrpaths) > 1:
        data = create_diff_table(data, out_path, comparison)
        data = fisher_test_cci(data, 'LRScore', out_path, comparison)
        data = mannwhitneyu_test_cci(data, 'LRScore', out_path, comparison)
        data = filtered_graphs(data, out_path)

    print("Calculating CCI Ranking")
    data = ranking(data, out_path, sel_columns=sel_columns, slot="graphs")
    print("Calculating GCI Ranking")
    data = ranking(data, out_path, sel_columns=sel_columns, slot="graphs_ggi")
    print("Network Analysis Done")

    with open(os.path.join(out_path, "LR_data.pkl"), "wb") as f:
        pickle.dump(data, f)

    return(data)

