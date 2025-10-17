import os
import pandas as pd
from tqdm import tqdm
from anndata import AnnData
from .Comparative_condition import *
from .Single_Condition import *
from .utils import *


def analise_LR(input,
                genes=None, 
                tf_genes=None, 
                out_path=None, 
                sep=',', 
                threshold=0, 
                colors=None, 
                out_file=None, 
                output_fmt="html_document", 
                sel_columns=['source','target','gene_A','gene_B','type_gene_A','type_gene_B','MeanLR'], 
                org='hsa', comparison=None, filtered_net=False, filename=None):
    
    """
    Core engine to generate report. Here we perform all the computation related to pyCrossTalkeR

    Parameters
    ----------
    input :
        Named vector with the lrpaths of each output or an AnnData object    
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
    filename :
        filename prefix for output files (to be provided if not already present in AnnData object)
    
    Returns
    -------
    Rmarkdown report all objects from each step
    
    """
    
    annData = read_lr_single_condition(input,
                                    sel_columns, 
                                    out_path, 
                                    sep, 
                                    colors)

    if isinstance(input, dict) or (isinstance(input, AnnData) and 'details' not in annData.uns['pycrosstalker']):
        annData.uns['pycrosstalker']['details'] = {}
        annData.uns['pycrosstalker']['details']['filename'] = filename

    print("Create a Differential Table")
    if len(annData.uns['pycrosstalker']['path']) > 1:
        annData = create_diff_table(annData, out_path, comparison)
        annData = fisher_test_cci(annData, 'LRScore', out_path, comparison)
        annData = mannwhitneyu_test_cci(annData, 'LRScore', out_path, comparison)
        annData = filtered_graphs(annData, out_path)

    print("Calculating CCI Ranking")
    annData = ranking(annData, out_path, sel_columns=sel_columns, slot="graphs")
    print("Calculating GCI Ranking")
    annData = ranking(annData, out_path, sel_columns=sel_columns, slot="graphs_ggi")
    print("Network Analysis Done")

    print("Generating h5ad file with Analysed Results")
    annData.write(os.path.join(out_path, annData.uns['pycrosstalker']['details']['filename'] + "_analysed.h5ad"))

    return(annData)

