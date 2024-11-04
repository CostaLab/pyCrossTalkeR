#'
#'Run and Generate all LR Downstream analysis
#'
#'This function loads the single conditions LR outputs and return the LR network
#`based analysis.
#'It assumes that the table present the following columns Ligand,
#` Ligand.Cluster,Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'@slot graphs All Cell Cell Interaction Networks
#'@slot tables All tables from single condition
#'@slot max_iter  Max meanLR from all
#'@slot max_nodes All Celltype in the experiment
#'@slot coords  Cell Cell Interaction Plots
#'@slot colors  Cell type colors
#'@slot rankings Ranking of cells and Genes
#'@slot loadings  CCI values to remove multiple times genes
#'@slot pca  PCA results
#'@slot annot  Annotation Results


import attr

@attr.s
class LRObj:
    """
    This function loads the single conditions LR outputs and return the LR network It assumes that the table present the following columns Ligand, measure

    Attributes
    ----------
    graphs :
        All Cell Cell Interaction Networks
    tables :
        All tables from single condition
    max_iter :
        Max meanLR from all
    max_nodes :
        All Celltype in the experiment
    coords :
        Cell Cell Interaction Plots
    colors :
        Cell type colors
    rankings :
        Ranking of cells and Genes
    loadings :
        CCI values to remove multiple times genes
    pca :
        PCA results
    annot :
        Annotation Results
    
    """

    graphs = attr.ib(default=attr.Factory(dict))
    graphs_ggi = attr.ib(default=attr.Factory(dict))
    tables = attr.ib(default=attr.Factory(dict))
    max_iter = attr.ib(default=0)
    max_nodes = attr.ib(default=0)
    coords = attr.ib(default=None)
    colors = attr.ib(default=None)
    rankings = attr.ib(default=attr.Factory(list))
    loadings = attr.ib(default=attr.Factory(list))
    pca = attr.ib(default=attr.Factory(list))
    annot = attr.ib(default=attr.Factory(list))
    stats = attr.ib(default=attr.Factory(list))
