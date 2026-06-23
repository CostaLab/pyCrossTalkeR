import pandas as pd
import numpy as np
import networkx as nx
import igraph as ig
import leidenalg
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
import seaborn as sns
import plotly.graph_objects as go
from adjustText import adjust_text
from gprofiler import GProfiler
from sankeyflow import Sankey

def plot_cci(graph, colors, plt_name, coords, pg, emax=None, leg=False, low=25, high=75, ignore_alpha=False, log=False, efactor=8, vfactor=12, vnames=True, figsize=None, scale_factor=2, node_size=2, font_size=10,return_figure=False):
    """
    This function does a CCI plot

    Parameters
    ----------
    graph :
        Paths of single condition LR data
    colors :
        Cell type (Cluster) Colors
    plt_name :
        Plot Name (Title)
    coords :
        object coordinates
    emax :
        Max MeanLR across the all inputs, if its not defined, the method going to consider the max find within a sample
    leg :
        Set color legend
    low :
        Lower threshold: This parameter low and high defines the edges
    high :
        Higher threshould which will be filtered. Edges within the interval [low\,high] are filtered.
    ignore_alpha :
        Not include transparency on the plot edges
    log :
        Logscale the interactions
    efactor :
        Edge scale factor
    vfactor :
        Certex scale factor
    vnames :
        Remove vertex labels
    pg :
        Pagerank values
    figsize:
        Set matplotlib figsize
    return_figure:
        Option for return matplotlib figure axes

    Returns
    -------
    Python default plot
    
    """
    graph = nx.from_pandas_edgelist(graph,
                                    source='source',
                                    target='target',
                                    edge_attr=True,
                                    create_using=nx.DiGraph())

    # Check Maximal Weight
    if emax is None:
        emax = 0
        for _, _, d in graph.edges(data=True):
            weight = d.get('weight', 0)
            if isinstance(weight, (int, float)):
                emax = max(emax, abs(weight))

    # Create colormap
    colors_list = sns.color_palette("coolwarm", 201)  # Adjust to match the R colormap
    col_pallet_colors = [colors_list[i] for i in range(201)]
    col_pallet_colors[10] = '#B8b9ba'  # Expand the palette range

    # Scale coordinates
    coords_array = np.array(list(coords.values()))

    if coords_array.shape[0] != 1:
        coords_mean = (np.mean(coords_array[:, 0]), np.mean(coords_array[:, 1]))
        coords_std = (np.std(coords_array[:, 0]), np.std(coords_array[:, 1]))
        coords_scale = {node: tuple((coord - mean) / std for coord, mean, std in zip(coords[node], coords_mean, coords_std)) for node in coords}
    else:
        coords_scale = coords

    coords_scale = {key: (coords_scale[key][0] * scale_factor, coords_scale[key][1] * scale_factor) for key in coords_scale}

    # Calculate edge colors and alpha
    edge_colors = []
    alpha = []

    for u, v, d in graph.edges(data=True):
        weight = d.get('weight', 0)
        we = np.round(np.interp(weight, [-emax, emax], [1, 200]))
        edge_colors.append(col_pallet_colors[int(we)])
        alpha_cond = low < d.get('inter', 0) < high and not np.isnan(d.get('inter', 0) < high)
        alpha.append(0 if alpha_cond else d.get('inter', 0) < high)

    # Set edge attributes
    for u, v, d in graph.edges(data=True):
        d['color'] = [(c[0], c[1], c[2], a) for c, a in zip(edge_colors, alpha)]
        if log:
            d['width'] = np.log2(1 + d.get('inter', 0)) * efactor if d.get('inter', 0) != 0 else 0
        else:
            d['width'] = d.get('inter', 0) * efactor if d.get('inter', 0) != 0 else 0
        d['arrow_size'] = 0.4
        d['arrow_width'] = d['width'] + 0.8
        d['loop_angle'] = np.nan

    node_colors = [str(colors.get(node)) for node in graph.nodes()]
    node_sizes = [size*1000*node_size for size in pg]

    # Plot the graph
    if figsize is None:
        figsize = (6,6)
    fig, ax = plt.subplots(figsize=figsize)
    
    nx.draw(graph, pos=coords_scale, edge_color=edge_colors, node_color=node_colors, node_size=node_sizes,
            width=[d['width'] for _, _, d in graph.edges(data=True)],
            arrows=True, arrowsize=30, arrowstyle='-|>',
            connectionstyle='arc3,rad=0.3', with_labels=vnames, font_size=font_size)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)

    # Node Pagerank legend
    if pg is not None:
        min_pg, max_pg = min(pg), max(pg)
        legend1 = ax.legend(loc='lower left', title="Pagerank",
                handles=[plt.Line2D([], [], linestyle='', marker='o', markersize=v / vfactor, markerfacecolor='black', markeredgecolor='none') for v in [min_pg, (min_pg + max_pg) / 2, max_pg]],
                labels=[round(min_pg, 2), round((min_pg + max_pg) / 2, 2), round(max_pg, 2)], bbox_to_anchor=(0.95, 0.3))

    # Thickness legend
    non_zero_inter_edges = [d['inter'] for _, _, d in graph.edges(data=True) if d.get('inter', 0) != 0]
    if non_zero_inter_edges:
        e_wid_sp = [round(min(non_zero_inter_edges), 2), round(min(non_zero_inter_edges) + (emax / 2), 2), round(emax, 2)]
        legend2 = ax.legend(e_wid_sp, title='Percentage of \nthe interactions', title_fontsize='small', loc='upper left', bbox_to_anchor=(0.95, 0.7))

    ax.add_artist(legend1)
    ax.set_title(plt_name)
    # Show the plot
    plt.tight_layout()
    plt.show()
    if return_figure:
        return (fig,ax)


def plot_pca_LR_comparative(lrobj_tblPCA, pca_table, dims=(1, 2), ret=False, ggi=True, include_tf=False, gene_types="all"):
    """
    This function is a proxy to the PCA plot in comparative conditions

    Parameters
    ----------
    lrobj_tblPCA :
        LRobject table with all data
    pca_table :
        table entry
    dims :
        PCA dims
    ret :
        return plot
    ggi :
        GGI mode
    include_tf :
        intracellular option
    gene_types :
        filter option of genes
    
    Returns
    -------
    Python default plot
    
    """
    
    pca_plot = {}
    # Extract PCA results and create a DataFrame
    pca_result = lrobj_tblPCA['pca'][pca_table]
    pca_df = lrobj_tblPCA['rankings'][pca_table]
    pca_df[['PC1', 'PC2']] = pca_result # Adjust dims to be zero-indexed
    pca_df = pca_df.set_index("nodes")

    if ggi:
        # Filter for LR or TF
        if gene_types == "LR":
            result_split_names = [name for name in pca_df.index if "|R" in name or "|L" in name]
        elif gene_types == "TF":
            result_split_names = [name for name in pca_df.index if "|TF" in name]
        else:
            result_split_names = pca_df.index.tolist()

        pca_df = pca_df.loc[result_split_names]

        # Mapping Table
        if include_tf:
            map_df = pd.DataFrame(pca_df.index, columns=["gene"])
            map_df["mapping"] = map_df["gene"].apply(lambda gene: "Receptor" if "|R" in gene else ("Ligand" if "|L" in gene else "Transcription Factor"))
            color_groups = ["#f8756b", "#00b835", "#619cff"]
        else:
            l_mapping = lrobj_tblPCA["tables"][pca_table.replace('_ggi', '')][["ligpair", "type_gene_A"]].rename(columns={"ligpair": "gene", "type_gene_A": "mapping"}).drop_duplicates()
            r_mapping = lrobj_tblPCA["tables"][pca_table.replace('_ggi', '')][["recpair", "type_gene_B"]].rename(columns={"recpair": "gene", "type_gene_B": "mapping"}).drop_duplicates()
            map_df = pd.concat([l_mapping, r_mapping]).drop_duplicates().reset_index(drop=True)
            map_df = map_df[map_df["gene"].isin(pca_df.index)]
            color_groups = ["#f8756b", "#00b835"]

        # Merge mapping info with PCA data
        pca_df = pca_df.merge(map_df, left_index=True, right_on="gene")

        #Threshold

        sdev_x = pca_df['PC1'].std()
        sdev_y = pca_df['PC2'].std()
        ver_zx = np.abs(pca_df['PC1']) >= (4 * sdev_x)
        ver_zy = np.abs(pca_df['PC2']) >= (4 * sdev_y)

        # Plotting
        x_max = max(abs(pca_df['PC1']))
        y_max = max(abs(pca_df['PC2']))

        pca_df['PC1'] = -pca_df['PC1']
        plt.figure(figsize=(10, 7))
        sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=20, hue='mapping', palette=color_groups)

        # Adjust text labels to avoid overlap
        texts = []
        for i, gene in enumerate(pca_df['gene']):
            if ver_zx.iloc[i] or ver_zy.iloc[i]: 
                texts.append(plt.text(pca_df.loc[pca_df['gene'] == gene, 'PC1'].values[0], pca_df.loc[pca_df['gene'] == gene, 'PC2'].values[0], gene, fontsize=8,  bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')))
        
        adjust_text(texts,
                    force_text=(1.0, 2.0))

        plt.xlim(-x_max, x_max)
        plt.ylim(-y_max, y_max)
        plt.xlabel(f'PC{dims[0]}')
        plt.ylabel(f'PC{dims[1]}')
        plt.title(pca_table, y=1.08)
        plt.legend(title='Gene Type')
        plt.grid(True)

        plt.axhline(0, linestyle='--', color='gray')
        plt.axvline(0, linestyle='--', color='gray')
        
        plt.show()

        pca_plot[pca_table] = plt
        
    else:
        # No GGI
        x_max = max(abs(pca_df['PC1']))
        y_max = max(abs(pca_df['PC2']))

        pca_df['PC1'] = -pca_df['PC1']
        plt.figure(figsize=(10, 7))
        sns.scatterplot(x='PC1', y='PC2', data=pca_df)

        # Adjust text labels to avoid overlap
        texts = []
        for i, gene in enumerate(pca_df.index):
            texts.append(plt.text(pca_df.loc[gene, 'PC1'], pca_df.loc[gene, 'PC2'], gene, fontsize=8))
        
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

        plt.xlim(-x_max, x_max)
        plt.ylim(-y_max, y_max)
        plt.xlabel(f'PC{dims[0]}')
        plt.ylabel(f'PC{dims[1]}')
        plt.title(pca_table)
        plt.grid(True)

        # Set x and y axis intervals
        plt.xticks(np.arange(-np.ceil(x_max), np.ceil(x_max) + 1))
        plt.yticks(np.arange(-np.ceil(y_max), np.ceil(y_max) + 1))

        plt.axhline(0, linestyle='--', color='gray')
        plt.axvline(0, linestyle='--', color='gray')
        plt.show()

        pca_plot[pca_table] = plt

    if ret:
        return pca_plot
    

def plot_bar_rankings(annData, table_name, ranking, type = None, filter_sign = None, mode = "cci", top_num = 10):
    """
    This function generates the barplot for a given network ranking on the CGI level. Further, the genes can be filtered by selected gene types to filter the plot.

    Parameters
    ----------
    annData :
        AnnData object with all data

    table_name :
        name of the ranking table

    ranking :
        name of the network ranking to use

    type :
        gene type (L,R,TF, LR/RL, RTF/TFR, LTF/TFL)

    filter_sign :
        show all (NULL), only positive (pos), or only negativ (neg) results
    
    Returns
    -------
    Python default plot

    """

    if '_x_' in table_name:
        rankings_table = annData.uns['pycrosstalker']['results']['rankings'][table_name]

        if type is not None:
            if len(type) == 1:
                rankings_table = rankings_table[rankings_table['nodes'].str.contains(r'\|' + type)]
            elif len(type) == 2:
                if type == 'TF':
                    rankings_table = rankings_table[rankings_table['nodes'].str.contains(r'\|' + type)]
                else:
                    rankings_table = rankings_table[rankings_table['nodes'].str.contains(r'\|R|\|L')]
            elif len(type) == 3:
                if type in ['RTF', 'TFR']:
                    rankings_table = rankings_table[rankings_table['nodes'].str.contains(r'\|R|\|TF')]
                elif type in ['LTF', 'TFL']:
                    rankings_table = rankings_table[rankings_table['nodes'].str.contains(r'\|L|\|TF')]

        rankings_table = rankings_table.sort_values(by=ranking)
        
        if mode == 'cgi':
            rankings_table = pd.concat([rankings_table.head(top_num), rankings_table.tail(top_num)])
            # rankings_table.loc[rankings_table[ranking].abs().nlargest(20).index]
        else:
            pass

        if filter_sign == 'pos':
            rankings_table = rankings_table[rankings_table['ranking'] > 0]
        elif filter_sign == 'neg':
            rankings_table = rankings_table[rankings_table['ranking'] < 0]


        if rankings_table.empty:
            return "No entries with provided Filters."

        rankings_table['signal'] = ['negative' if x < 0 else 'positive' for x in rankings_table[ranking]]

        custom_palette = {'positive': '#FF6E00', 'negative':'#00FFFF'}  # Orange and Blue

        # Plot
        plt.figure(figsize=(8, 6))
        sns.barplot(x=ranking, y='nodes', data=rankings_table, hue='signal', dodge=False, palette=custom_palette)
        plt.title(f"Ranking for {table_name}")
        plt.xlabel(ranking)
        plt.ylabel('Nodes')

        # Set x-axis tick intervals
        max_val = rankings_table[ranking].max()
        min_val = rankings_table[ranking].min()
        ticks = np.linspace(min_val, max_val, num=5)  # Adjust 'num' for more/less intervals
        plt.xticks(ticks, [f'{tick:.2f}' for tick in ticks])

        # Invert y-axis to have highest values at the top
        plt.gca().invert_yaxis()

        # Show the legend only once
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles, labels, loc='lower right')

        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.gca().set_axisbelow(True)
        # Adjust layout and show plot
        plt.tight_layout()
        plt.show()


def plot_sankey(lrobj_tbl, target = None, ligand_cluster = None, receptor_cluster = None, plt_name = None, threshold = 50, tfflag = True):
    """
    This function selected genes sankey plot

     Parameters
    ----------
    lrobj_tbl :
        LRobject table with all data

    target :
        gene

    ligand_cluster :
        Ligand Clusters

    receptor_cluster :
        Receptor Clusters

    plt_name :
        plot title

    threshold :
        top_n n value
    
    Returns
    -------
    Python default plot

    """

    lrobj_tbl = lrobj_tbl[(lrobj_tbl['type_gene_A'] == "Ligand") & (lrobj_tbl['type_gene_B'] == "Receptor")]

    if target is not None:
        if len(target.split('|')) > 1:
            target_type = str(target.split('|')[1])
            if target_type == 'R':
                if lrobj_tbl['gene_B'].str.contains('\\|').any():
                    pass
                else:
                    target = target.split('|')[0]
                data = lrobj_tbl[lrobj_tbl['gene_B'] == target]
            elif target_type == 'L':
                if lrobj_tbl['gene_A'].str.contains('\\|').any():
                    pass
                else:
                    target = target.split('|')[0]
                data = lrobj_tbl[lrobj_tbl['gene_A'] == target]
        else:
            data = lrobj_tbl[lrobj_tbl['allpair'].str.contains(target)]
    else:
        data = lrobj_tbl

    
    if ligand_cluster is not None:
        data = data[data['source'].isin(ligand_cluster)]
    
    if receptor_cluster is not None:
        data = data[data['target'].isin(receptor_cluster)]

    color_palette = ['#00BFC4', '#FF3E3E']

    
    if len(data) >= 1:
        cat_cols = ['source', 'gene_A', 'gene_B', 'target']
        value_cols = 'LRScore'
        data = data.loc[data['LRScore'].abs().nlargest(min(len(data), threshold)).index]
        title = plt_name

        gen_sankey(data, cat_cols, value_cols, title)
    
    else:
        print(f"Gene->{target} Not Found")
    

def gen_sankey2(df, cat_cols=[], value_cols='', title='Sankey Diagram'):
    """
    Helper function to the function plot_sankey()

     Parameters
    ----------
    df :
        Dataframe

    cat_cols :
        Columns interested in the sankey plot

    value_cols :
        Sankey plot generated using connections based on this value_cols

    title :
        Title of Sankey plot
    
    Returns
    -------
    Nothing (plots Sankey plot)

    """

    df['source'] += 'S'
    df['target'] += 'T'
    
    labelList = []
    for catCol in cat_cols:
        labelListTemp =  list((df[catCol].values))
        labelList = labelList + labelListTemp    
        
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            tempDf.columns = ['source','target','count']
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
        # sourceTargetDf = sourceTargetDf.groupby(['source','target']).agg({'count':'sum'}).reset_index()
        
    sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
    sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))

    for i, label in enumerate(labelList):
        if label[-1:] == 'S' or label[-1:] == 'T':
            labelList[i] = labelList[i][:-1]
        
    norm = mcolors.Normalize(vmin=min(sourceTargetDf['count']), vmax=max(sourceTargetDf['count']))
    colormap = cm.get_cmap('RdBu_r')
    link_colors = [mcolors.to_hex(colormap(norm(value))) for value in sourceTargetDf['count']]
    
    fig =  go.Figure(data = [go.Sankey(
        node = dict(
          pad = 0,
          thickness = 20,
          line = dict(
            color = "black",
            width = 0.5
          ),
          label = labelList,
          color = "white"
        ),
        link = dict(
          source = sourceTargetDf['sourceID'],
          target = sourceTargetDf['targetID'],
          value = [abs(i) for i in sourceTargetDf['count']],
          color = link_colors
        )    
    )])
    
    colorbar_trace = go.Scatter(
        x=[None], y=[None], mode='markers',
        marker=dict(
            colorscale='RdBu_r',
            cmin=min(sourceTargetDf['count']),
            cmax=max(sourceTargetDf['count']),
            colorbar=dict(
                title="Value",
                thickness=15,
                len=0.5,
                x=1.05,
                xref="paper"
            )
        ),
        hoverinfo='none'
    )

    fig.add_trace(colorbar_trace)

    fig.add_annotation(x=0, y=1.05, yref="paper", text="Source", showarrow=False, font=dict(size=10))
    fig.add_annotation(x=0.33, y=1.05, yref="paper", text="Ligand", showarrow=False, font=dict(size=10))
    fig.add_annotation(x=0.66, y=1.05, yref="paper", text="Receptor", showarrow=False, font=dict(size=10))
    fig.add_annotation(x=1, y=1.05, yref="paper", text="Target", showarrow=False, font=dict(size=10))

    fig.update_layout(
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0,1]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0,1]),
        plot_bgcolor='white',
        autosize=True,
        width = None,
        height = 600,
        title = title,
        font = dict(size=10)
        )
    
    fig.show(config={"responsive": True})

def gen_sankey(df, cat_cols=[], value_cols='', title='Sankey Diagram'):
    """
    Helper function to the function plot_sankey()

     Parameters
    ----------
    df :
        Dataframe

    cat_cols :
        Columns interested in the sankey plot

    value_cols :
        Sankey plot generated using connections based on this value_cols

    title :
        Title of Sankey plot
    
    Returns
    -------
    Nothing (plots Sankey plot)

    """

    # df['source'] += 'S'
    df['target'] += ' '
    
    labelList = []
    for catCol in cat_cols:
        labelListTemp =  list((df[catCol].values))
        labelList = labelList + labelListTemp    
        
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            tempDf.columns = ['source','target','count']
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
    
    vmin = sourceTargetDf['count'].min()
    vmax = sourceTargetDf['count'].max()
    limit = max(abs(vmin), abs(vmax))
    vcenter = 0
    norm = mcolors.TwoSlopeNorm(vmin=-limit, vcenter=vcenter, vmax=limit)
    # norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap('RdBu_r')
    sourceTargetDf['hex_color'] = sourceTargetDf['count'].apply(lambda x: mcolors.to_hex(cmap(norm(x))))
    
    flows = []
    for i, row in sourceTargetDf.iterrows():
        flows.append((row['source'], row['target'], 1, {'color': row['hex_color']}))

    nodes = Sankey.infer_nodes(flows)
    nodes_new = []
    for level in nodes:
        level_new = []
        for node in level:
            node_new = node + [{'color' : 'black',
                                'label_pos':'center', 'label_opts': dict(fontsize=10, bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))}]
            level_new.append(node_new)
        nodes_new.append(level_new)

    fig, ax = plt.subplots(figsize=(15, 10))
    s = Sankey(flows=flows,
               nodes=nodes_new,
               flow_color_mode_alpha=0.3,
               node_opts=dict(label_format='{label}'),
    )
    s.draw(ax=ax)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.01, shrink=0.5)
    cbar.set_label(value_cols, fontsize=12)

    ax.text(x=-0.05, y=1.02, s="Source", fontsize=10)
    ax.text(x=0.95, y=1.02, s="Ligand", fontsize=10)
    ax.text(x=1.95, y=1.02, s="Receptor", fontsize=10)
    ax.text(x=2.95, y=1.02, s="Target", fontsize=10)

    ax.set_title(title)
    plt.tight_layout()
    plt.show()

def gene_annotation(gene_list_to_profile,
                    num_gos: int = 15,
                    figsize=(10,6),
                    title: str = None,
                    font_size: int = 10,
                    organism: str = 'hsapiens',
                    dpi: int = 100,
                    s: int = 100,
                    color: str = 'tab:blue'):
    """
    Perform Gene Ontology (GO) enrichment analysis and create a scatterplot of enriched terms.

    Parameters:
    ----------
    num_gos: int, optional
        Number of GO terms to plot. Default is 5.
    figsize: tuple, optional
        figsize. Default is (6,6).
    title: str
        Title of the plot.
    font_size: int, optional
        Font size for labels. Default is 10.
    0rganism: str, optional
        The organism for GO analysis. Default is 'hsapiens'.
    dpi: int, optional
        Dots per inch for the saved plot image. Default is 100.
    s: int, optional
        Marker size for the scatterplot. Default is 200.
    color: str, optional
        Color of the scatterplot markers. Default is 'tab:blue'.

    Returns:
    --------
    None
        Plots the scatterplot of enriched GO terms.
    """

    
    gp = GProfiler(return_dataframe=True)
    if gene_list_to_profile:
        gprofiler_results = gp.profile(organism = organism,
                                       query = gene_list_to_profile)
    else:
        return "Genes list is empty!"
    
    
    if(gprofiler_results.shape[0] == 0):
        return "Not enough information!"

    
    if(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]

  
    selected_gps = gprofiler_results.head(num_gos)[['name', 'p_value']]
    
    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    plt.figure(figsize = figsize, dpi = dpi)
    # plt.style.use('default')
    sns.scatterplot(data = selected_gps, x = "nlog10", y = "name", s = s, color = color)

    plt.title(title, fontsize = font_size)

    plt.xticks(size = font_size)
    plt.yticks(size = font_size)

    plt.ylabel("GO Terms", size = font_size)
    plt.xlabel("-$log_{10}$ (P-value)", size = font_size)

    plt.tight_layout()
    plt.show()


def plot_volcane(df, method, p_threshold=0.05, fc_threshold=1, figsize=(8, 6), annot=True, title=None):
    """
    This function generates a Volcano plot

    Parameters
    ----------
    df :
        Dataframe

    Returns
    -------
    Python default volcano plot

    """
    np.random.seed(42)
    data = df
    data['neg_log10_p_value'] = -np.log10(df['p_value'])
    if method == "fisher":
        attr = "lodds"
    elif method == "mannwhitneyu":
        attr = "lfc"

    data["color"] = "gray"
    data.loc[(data[attr] > fc_threshold) & (data["p_value"] < p_threshold), "color"] = "red"
    data.loc[(data[attr] < -fc_threshold) & (data["p_value"] < p_threshold), "color"] = "red"
    data.loc[(data[attr] > -fc_threshold) & (data[attr] < fc_threshold) & (data["p_value"] < p_threshold), "color"] = "blue"
    data.loc[(data[attr] < -fc_threshold) & (data["p_value"] > p_threshold), "color"] = "green"
    data.loc[(data[attr] > fc_threshold) & (data["p_value"] > p_threshold), "color"] = "green"

    # Plot
    plt.figure(figsize=figsize)
    sns.scatterplot(x=attr, y="neg_log10_p_value", hue="color", palette={"gray": "gray", "red": "red", "blue": "blue", "green": "green"}, data=data, edgecolor=None, alpha=0.7)

    # Add significance threshold lines
    plt.axhline(-np.log10(p_threshold), linestyle="--", color="black", linewidth=1)
    plt.axvline(fc_threshold, linestyle="--", color="black", linewidth=1)
    plt.axvline(-fc_threshold, linestyle="--", color="black", linewidth=1)

    if annot:
        for i, row in data.iterrows():
            if row['color'] == 'red':
                plt.text(row[attr], row["neg_log10_p_value"], row["cellpair"], fontsize=8, ha='right')

    x_limit = max(abs(data[attr].min()), abs(data[attr].max()))
    plt.xlim(-x_limit-1, x_limit+1)

    plt.xlabel(r"Log$_{2}$ Fold Change")
    plt.ylabel(r"-Log$_{10}$(p-value)")
    plt.title(title)
    plt.legend([],[], frameon=False)
    plt.show()


def plot_clustermap(data, title, annot=True, return_figure=False):
    """
    This function generates a Clustermap plot

    Parameters
    ----------
    data :
        Dataframe
    title : str
        Title of the plot
    annot : bool
        Whether to annotate the heatmap with values
    return_figure:
        Option for return matplotlib figure axes
    Returns
    -------
    
    Python Cluster map

    """
    pivot_table = data.groupby(["source", "target"])["LRScore"].sum().unstack().fillna(0)
    xlabel, ylabel = "Target", "Source"
    g = sns.clustermap(
        pivot_table,
        figsize=(9, 7),
        annot=annot,
        linewidths=0.5,
        method="ward",
        metric="euclidean",
        dendrogram_ratio=(0.2, 0.2),
        cbar_pos=(0.02, 0.8, 0.03, 0.15)
    )
    g.ax_heatmap.set_xlabel(xlabel)
    g.ax_heatmap.set_ylabel(ylabel)
    plt.title(title, fontsize=14)
    plt.show()


def plot_graph_clustermap(graph, weight="LRScore", title="Ligand-Receptor Heatmap", annot=True):
    """
    This function generates the graph adjacency matrix Heatmap

    Parameters
    ----------
    data :
        Dataframe
    weight : str
        The weight attribute to use for the adjacency matrix, default is "LRScore"
    title : str
        Title of the plot
    annot : bool
        Whether to annotate the heatmap with values

    Returns
    -------
    Python Cluster map

    """
    graph = nx.from_pandas_edgelist(graph,
                                    source='source',
                                    target='target',
                                    edge_attr=True,
                                    create_using=nx.DiGraph())

    nodes = list(graph.nodes)
    adj_matrix = nx.to_pandas_adjacency(graph, nodelist=nodes, weight=weight).fillna(0).astype(float)
    max_val = np.abs(adj_matrix.values).max()

    g = sns.clustermap(
        adj_matrix,
        figsize=(9, 7),
        cmap="RdBu_r",
        linewidths=0.5,
        center=0,
        vmin=-max_val,
        vmax=max_val,
        annot=annot,
        method="ward",
        metric="euclidean",
        dendrogram_ratio=(0.2, 0.2),
        cbar_pos=(0.02, 0.8, 0.03, 0.15)
    )

    g.ax_heatmap.set_xlabel("Receptor Cluster")
    g.ax_heatmap.set_ylabel("Ligand Cluster")
    plt.title(title, fontsize=14)

    plt.show()

def parse_node(node_string):
    """
    This function parses a node string into its celltype, gene and class

    Parameters
    ----------
    node_string :
        Node label, e.g. "ImmatureNeutrophils/S100a8|L"

    Returns
    -------
    Dict with keys 'celltype', 'gene' and 'class'

    """

    celltype, rest = node_string.split("/", 1)

    gene, gene_class = rest.rsplit("|", 1)

    return {
        "celltype": celltype,
        "gene": gene,
        "class": gene_class
    }

def _prune_orphan_nodes(plot_nodes, plot_edges):
    """
    This function removes nodes not touched by any edge and rebuilds the layer index

    Parameters
    ----------
    plot_nodes :
        Mapping of node name to node attributes

    plot_edges :
        List of edges, each with 'source' and 'target' keys

    Returns
    -------
    Tuple of pruned plot_nodes and rebuilt step->nodes layer index

    """

    used_nodes = set()

    for edge in plot_edges:
        used_nodes.add(edge["source"])
        used_nodes.add(edge["target"])

    plot_nodes = {
        node: attrs
        for node, attrs in plot_nodes.items()
        if node in used_nodes
    }

    nodes_by_layer = defaultdict(list)

    for node in plot_nodes:
        step = int(node.split("__step")[-1])
        nodes_by_layer[step].append(node)

    return plot_nodes, nodes_by_layer


def _as_list(value):
    """
    This function normalizes a single value or a list/tuple/set into a plain list

    Parameters
    ----------
    value :
        Single value, list/tuple/set of values, or None

    Returns
    -------
    A plain list, or None if value is None

    """

    if value is None:
        return None

    if isinstance(value, (list, tuple, set)):
        return list(value)

    return [value]


def plot_intracellular_signaling_network(
    lrobj_tbl,
    mode="LRTFL",
    gene=None,
    level=None,
    celltype=None,
    ligand_celltype=None,
    receptor_celltype=None,
    font_size=10,
    node_spacing=2,
    celltype_spacing=6,
    step_spacing=4,
    edge_width=4,
    arrow_size=1.2,
    top_n_per_step=None,
    top_n_lr=None,
    top_n_rtfl=None,
    position="any"
):
    """
    This function plots multi-step ligand-receptor-transcription factor
    signaling cascades as a layered network, where each node is a
    gene in a celltype (ligand, receptor or TF) and each edge is a
    directed interaction coloured by its LRScore. The paths traced
    depend on `mode`, and can be narrowed to genes, celltypes or the
    strongest interactions via the filtering arguments.

    Parameters
    ----------
    lrobj_tbl : pd.DataFrame
        Signaling network data on gene level

    mode : str
        Path topology: "RTFL" (R->TF->L), "LRTFL" (L->R->TF->L) or "RTFLR" (R->TF->L->R)

    gene : str or None
        Gene of interest to filter on

    level : str or None
        Restrict gene filtering to a molecule class: "L", "R" or "TF"

    celltype : str or None
        Restrict gene/level filtering to a specific celltype (any node in the path)

    ligand_celltype : str, list of str, or None
        Restrict to paths whose source ligand belongs to the given celltype(s).
        Only valid in "LRTFL" mode

    receptor_celltype : str, list of str, or None
        Restrict to paths whose receptor (and its intracellular TF and output
        ligand) belongs to the given celltype(s)

    position : str
        Where the gene+level match must occur: "any", "first" or "last"

    font_size : int
        Node label font size

    node_spacing : float
        Vertical spacing between nodes of the same celltype within a layer

    celltype_spacing : float
        Vertical spacing between different celltypes within a layer

    step_spacing : float
        Horizontal spacing between signaling steps

    edge_width : float
        Line width for edges

    arrow_size : float
        Size of arrowheads

    top_n_per_step : int or None
        Keep only the strongest N edges per step (by |weight|), pruning orphan nodes

    top_n_lr : int or None
        Cap on the strongest ligand->receptor edges (anchor stage in "LRTFL",
        last stage in "RTFLR", ignored in "RTFL")

    top_n_rtfl : int or None
        Cap on the strongest receptor->TF and TF->ligand edges, counted per
        receptor celltype

    Returns
    -------
    plotly.graph_objects.Figure

    """

    paths = []

    G = nx.MultiDiGraph()

    for _, row in lrobj_tbl.iterrows():

        source = row["ligpair"]
        target = row["recpair"]

        source_info = parse_node(source)
        target_info = parse_node(target)

        # Add source node
        G.add_node(
            source,
            celltype=source_info["celltype"],
            gene=source_info["gene"],
            molecule_class=source_info["class"]
        )

        # Add target node
        G.add_node(
            target,
            celltype=target_info["celltype"],
            gene=target_info["gene"],
            molecule_class=target_info["class"]
        )

        # Add directed edge
        G.add_edge(
            source,
            target,
            weight=row["LRScore"],
            interaction_type=row["interaction_type"],
            cellpair=row["cellpair"]
        )
    
    if mode == "RTFL":

        receptors = [
            n
            for n, attrs in G.nodes(data=True)
            if attrs["molecule_class"] == "R"
        ]

        for receptor in receptors:

            for _, tf, edge1 in G.out_edges(receptor, data=True):

                if edge1["interaction_type"] != "RTF":
                    continue

                for _, ligand, edge2 in G.out_edges(tf, data=True):

                    if edge2["interaction_type"] != "TFL":
                        continue

                    paths.append([receptor, tf, ligand])

    elif mode == "LRTFL":

        ligands = [
            n
            for n, attrs in G.nodes(data=True)
            if attrs["molecule_class"] == "L"
        ]

        for ligand in ligands:

            for _, receptor, edge1 in G.out_edges(ligand, data=True):

                if edge1["interaction_type"] != "LR":
                    continue

                for _, tf, edge2 in G.out_edges(receptor, data=True):

                    if edge2["interaction_type"] != "RTF":
                        continue

                    for _, target_ligand, edge3 in G.out_edges(tf, data=True):

                        if edge3["interaction_type"] != "TFL":
                            continue

                        paths.append([ligand, receptor, tf, target_ligand])

    elif mode == "RTFLR":

        receptors = [
            n
            for n, attrs in G.nodes(data=True)
            if attrs["molecule_class"] == "R"
        ]

        for receptor in receptors:

            for _, tf, edge1 in G.out_edges(receptor, data=True):

                if edge1["interaction_type"] != "RTF":
                    continue

                for _, ligand, edge2 in G.out_edges(tf, data=True):

                    if edge2["interaction_type"] != "TFL":
                        continue

                    for _, target_receptor, edge3 in G.out_edges(ligand, data=True):

                        if edge3["interaction_type"] != "LR":
                            continue

                        paths.append([receptor, tf, ligand, target_receptor])

    else:
        raise ValueError("mode must be 'RTFL', 'LRTFL' or 'RTFLR'")

    if position not in ("any", "first", "last"):
        raise ValueError("position must be 'any', 'first', or 'last'")

    if gene is not None or level is not None or celltype is not None:

        filtered = []

        for path in paths:

            keep = False
            last_index = len(path) - 1

            for idx, node in enumerate(path):

                attrs = G.nodes[node]

                if gene is not None and attrs["gene"] != gene:
                    continue

                if level is not None and attrs["molecule_class"] != level:
                    continue

                if celltype is not None and attrs["celltype"] != celltype:
                    continue

                if position == "first" and idx != 0:
                    continue

                if position == "last" and idx != last_index:
                    continue

                keep = True
                break

            if keep:
                filtered.append(path)

        paths = filtered

    ligand_celltype_list = _as_list(ligand_celltype)

    if ligand_celltype_list is not None:

        if mode != "LRTFL":
            raise ValueError(
                "ligand_celltype filters the incoming/source ligand "
                "node, which only exists in LRTFL mode (RTFL and RTFLR "
                "paths start at a receptor). Use mode='LRTFL', or filter "
                "the output ligand instead via celltype/position='last'."
            )

        paths = [
            path for path in paths
            if G.nodes[path[0]]["celltype"] in ligand_celltype_list
        ]

    receptor_celltype_list = _as_list(receptor_celltype)

    if receptor_celltype_list is not None:

        receptor_idx = 1 if mode == "LRTFL" else 0
        tf_idx = receptor_idx + 1
        target_ligand_idx = receptor_idx + 2

        filtered = []

        for path in paths:

            r_celltype = G.nodes[path[receptor_idx]]["celltype"]

            if r_celltype not in receptor_celltype_list:
                continue

            tf_celltype = G.nodes[path[tf_idx]]["celltype"]
            target_l_celltype = G.nodes[path[target_ligand_idx]]["celltype"]

            if tf_celltype != r_celltype or target_l_celltype != r_celltype:
                continue

            filtered.append(path)

        paths = filtered

    if len(paths) == 0:
        raise ValueError("No valid paths found after filtering.")

    if top_n_lr is not None or top_n_rtfl is not None:
        
        cascade = {
            "LRTFL": {
                "stages": [("global", top_n_lr),
                           ("per_celltype", top_n_rtfl),
                           ("per_celltype", top_n_rtfl)],
                "auth": True,
            },
            "RTFL": {
                "stages": [("per_celltype", top_n_rtfl),
                           ("per_celltype", top_n_rtfl)],
                "auth": False,
            },
            "RTFLR": {
                "stages": [("per_celltype", top_n_rtfl),
                           ("per_celltype", top_n_rtfl),
                           ("global", top_n_lr)],
                "auth": True,
            },
        }[mode]
        stages = cascade["stages"]

        def _edge_weight(u, v):
            return list(G.get_edge_data(u, v).values())[0]["weight"]

        def _top_edges(edge_keys, n):
            # edge_keys: iterable of (u, v) graph-node-name tuples.
            keys = set(edge_keys)
            if n is None:
                return keys
            ranked = sorted(
                keys, key=lambda uv: abs(_edge_weight(*uv)), reverse=True
            )
            return set(ranked[:n])

        def _top_edges_per_celltype(edge_keys, n):
            keys = set(edge_keys)
            if n is None:
                return keys
            buckets = defaultdict(set)
            for u, v in keys:
                buckets[G.nodes[u]["celltype"]].add((u, v))
            kept = set()
            for _, bucket in buckets.items():
                kept |= _top_edges(bucket, n)
            return kept

        def _select(k, scope, n):
            keys = ((q[k], q[k + 1]) for q in paths)
            if scope == "global":
                return _top_edges(keys, n)
            return _top_edges_per_celltype(keys, n)

        # --- Stage 0: the anchor ---
        scope0, n0 = stages[0]
        survivors = _select(0, scope0, n0)
        paths = [p for p in paths if (p[0], p[1]) in survivors]
        
        authoritative = cascade["auth"] and n0 is not None and len(stages) > 1

        reachable = {p[1] for p in paths} if authoritative else set()

        for k in range(1, len(stages)):

            scope_k, n_k = stages[k]
            survivors = _select(k, scope_k, n_k)

            if authoritative:
                best_per_source = {}
                for p in paths:
                    src_node, dst_node = p[k], p[k + 1]
                    if src_node not in reachable:
                        continue
                    cur = best_per_source.get(src_node)
                    if cur is None or abs(_edge_weight(src_node, dst_node)) > abs(_edge_weight(*cur)):
                        best_per_source[src_node] = (src_node, dst_node)
                guaranteed = set(best_per_source.values())
                survivors |= guaranteed
                reachable = {dst for (_, dst) in guaranteed}

            paths = [p for p in paths if (p[k], p[k + 1]) in survivors]

        if len(paths) == 0:
            raise ValueError(
                "No complete paths survived top-N filtering. Try raising "
                "top_n_lr / top_n_rtfl, or loosening other filters."
            )

    plot_nodes = {}
    plot_edges = []

    nodes_by_layer = defaultdict(list)

    for path in paths:

        for step, node in enumerate(path):

            plot_node = f"{node}__step{step}"

            if plot_node not in plot_nodes:

                attrs = G.nodes[node]

                plot_nodes[plot_node] = {
                    "celltype": attrs["celltype"],
                    "gene": attrs["gene"],
                    "molecule_class": attrs["molecule_class"],
                    "step": step
                }

                nodes_by_layer[step].append(plot_node)

        for step in range(len(path) - 1):

            source = path[step]
            target = path[step + 1]

            source_plot = f"{source}__step{step}"
            target_plot = f"{target}__step{step+1}"

            edge_data = list(G.get_edge_data(source, target).values())[0]

            plot_edges.append({
                "source": source_plot,
                "target": target_plot,
                "interaction_type": edge_data["interaction_type"],
                "weight": edge_data["weight"]
            })

    _seen_transitions = {}
    for edge in plot_edges:
        _seen_transitions[(edge["source"], edge["target"])] = edge
    plot_edges = list(_seen_transitions.values())

    if top_n_per_step is not None:

        step_edges = defaultdict(list)

        for edge in plot_edges:
            step = int(edge["source"].split("__step")[-1])
            step_edges[step].append(edge)

        filtered_edges = []

        for step, edges in step_edges.items():

            edges = sorted(edges, key=lambda x: abs(x["weight"]), reverse=True)
            filtered_edges.extend(edges[:top_n_per_step])

        plot_edges = filtered_edges

    # Drop incomplete signaling paths

    edge_lookup = {
        (edge["source"], edge["target"]): edge for edge in plot_edges
    }

    kept_transitions = set(edge_lookup.keys())
    final_transitions = set()

    for path in paths:

        transitions = [
            (f"{path[step]}__step{step}", f"{path[step + 1]}__step{step + 1}")
            for step in range(len(path) - 1)
        ]

        if all(t in kept_transitions for t in transitions):
            final_transitions.update(transitions)

    plot_edges = [edge_lookup[t] for t in final_transitions]

    if len(plot_edges) == 0:
        raise ValueError(
            "No complete signaling paths survived filtering -- every "
            "path lost at least one transition (e.g. its LR edge "
            "passed top_n_lr but its RTF/TFL edges didn't pass "
            "top_n_rtfl, or vice versa). Try raising top_n_lr / "
            "top_n_rtfl / top_n_per_step, or loosening other filters."
        )

    # Always prune dangling/orphaned nodes after any edge-level filtering,
    # regardless of which filter(s) above produced them.
    plot_nodes, nodes_by_layer = _prune_orphan_nodes(plot_nodes, plot_edges)

    if len(plot_nodes) == 0:
        raise ValueError("No nodes left after filtering -- thresholds are too strict.")

    # Node Positions
    # 1. compute each layer's internal y-positions and total height
    # 2. center every layer against the TRUE global max height
 
    layer_positions = {}
    layer_heights = {}

    for step, layer_nodes in nodes_by_layer.items():

        layer_nodes = sorted(
            layer_nodes,
            key=lambda n: (plot_nodes[n]["celltype"], plot_nodes[n]["gene"])
        )

        current_y = 0
        prev_celltype = None
        temp = []

        for node in layer_nodes:

            celltype_ = plot_nodes[node]["celltype"]

            if prev_celltype is not None:
                if celltype_ != prev_celltype:
                    current_y += celltype_spacing
                else:
                    current_y += node_spacing

            temp.append((node, current_y))
            prev_celltype = celltype_

        layer_positions[step] = temp
        layer_heights[step] = current_y

    global_height = max(layer_heights.values()) if layer_heights else 0

    pos = {}

    for step, temp in layer_positions.items():

        y_offset = (global_height - layer_heights[step]) / 2

        for node, y in temp:
            pos[node] = (step * step_spacing, y + y_offset)

    weights = np.array([e["weight"] for e in plot_edges])

    abs_max = max(np.max(np.abs(weights)), 1e-6)

    norm = mcolors.TwoSlopeNorm(vmin=-abs_max, vcenter=0, vmax=abs_max)

    cmap = cm.get_cmap("bwr")

    edge_traces = []
    arrow_annotations = []

    for edge in plot_edges:

        x0, y0 = pos[edge["source"]]
        x1, y1 = pos[edge["target"]]

        color = mcolors.to_hex(cmap(norm(edge["weight"])))

        edge_traces.append(
            go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(color=color, width=edge_width),
                hoverinfo="text",
                text=(
                    f"{edge['interaction_type']}<br>"
                    f"Score: {edge['weight']:.3f}"
                ),
                showlegend=False
            )
        )

        arrow_annotations.append(
            dict(
                x=x1, y=y1,
                ax=x0, ay=y0,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=3,
                arrowsize=arrow_size,
                arrowwidth=edge_width * 0.5,
                arrowcolor=color
            )
        )

    molecule_colors = {
        "L": "#2E8B57",
        "R": "#D62728",
        "TF": "#9467BD"
    }

    node_x = []
    node_y = []
    node_text = []
    node_hover = []
    node_color = []

    for node, attrs in plot_nodes.items():

        x, y = pos[node]

        node_x.append(x)
        node_y.append(y)

        node_text.append(
            f"{attrs['celltype']}<br>"
            f"{attrs['gene']} ({attrs['molecule_class']})"
        )

        node_hover.append(
            f"Celltype: {attrs['celltype']}<br>"
            f"Gene: {attrs['gene']}<br>"
            f"Class: {attrs['molecule_class']}"
        )

        node_color.append(molecule_colors.get(attrs["molecule_class"], "gray"))

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=node_text,
        textposition="top center",
        hovertext=node_hover,
        hoverinfo="text",
        marker=dict(
            size=10,
            color=node_color,
            line=dict(width=1, color="black")
        ),
        textfont=dict(size=font_size),
        showlegend=False
    )


    legend_names = {"L": "Ligand", "R": "Receptor", "TF": "Transcription Factor"}

    legend_traces = []

    for mol_type, color in molecule_colors.items():

        legend_traces.append(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=10, color=color),
                name=legend_names[mol_type],
                showlegend=True
            )
        )

    colorbar_trace = go.Scatter(
        x=[None],
        y=[None],
        mode="markers",
        marker=dict(
            size=0.1,
            color=[-abs_max, abs_max],
            colorscale="RdBu",
            cmin=-abs_max,
            cmax=abs_max,
            showscale=True,
            colorbar=dict(
                title="Score",
                thickness=15,
                len=0.35,
                y=0.2,
                x=1.02,
                outlinewidth=0.5,
                tickmode="array",
                tickvals=[-abs_max, -abs_max / 2, 0, abs_max / 2, abs_max],
                ticktext=[
                    f"{-abs_max:.2f}",
                    f"{-abs_max/2:.2f}",
                    "0",
                    f"{abs_max/2:.2f}",
                    f"{abs_max:.2f}"
                ]
            )
        ),
        hoverinfo="none",
        showlegend=False
    )

    fig = go.Figure(
        data=edge_traces + [node_trace] + legend_traces + [colorbar_trace]
    )

    fig.update_layout(
        width=2400,
        height=max(1200, int(global_height * 20)),
        plot_bgcolor="white",
        annotations=arrow_annotations,
        xaxis=dict(
            title="Signaling Step",
            showgrid=True,
            zeroline=False
        ),
        yaxis=dict(
            showgrid=False,
            showticklabels=False,
            zeroline=False
        ),
        legend=dict(
            title="Node Type",
            x=1.02,
            y=0.95,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="lightgray",
            borderwidth=1,
            font=dict(size=10)
        ),
        margin=dict(l=80, r=280, t=60, b=60)
    )

    return fig

def get_community_cmap(n):
    return plt.get_cmap("tab20", n)

def cci_community_layout(
    CCI_table,
    method="leiden",
    res=1.0,
    compression=0.2,
    seed=12,
    spring_k=0.1,
    spring_iterations_community=1000,
    spring_iterations_nodes=500,
    edge_width_scale=3.0,
    edge_alpha_intra=0.7,
    edge_alpha_inter=0.7,
    node_size=300,
    arrow_size=20,
    figsize=[10, 10],
):
    """
    Unified community layout for cell-cell communication networks.
    Parameters
    ----------
    method : str
        Community detection algorithm: "leiden" or "louvain".
    spring_iterations_nodes : int
        Only used when method="leiden" (per-community spring layout iterations).
    spring_k : float
        Spring layout k parameter. Leiden default: 0.1, Louvain default: 75.
    """
    # ---- Build similarity matrix ----
    W = CCI_table.pivot_table(
        index="from", columns="to", values="weight", fill_value=0.0
    )
    cells = W.index.to_numpy()

    deg = W.abs().sum(axis=1)
    deg_inv_sqrt = 1.0 / np.sqrt(deg.replace(0, np.nan))
    W_deg = W.mul(deg_inv_sqrt, axis=0).mul(deg_inv_sqrt, axis=1).fillna(0)
    K = 0.5 * (W_deg.T @ W_deg + W_deg @ W_deg.T)
    K_np = K.values

    # ---- Build undirected similarity graph ----
    G = nx.Graph()
    G.add_nodes_from(cells)
    for i in range(len(cells)):
        for j in np.argsort(K_np[i])[::-1]:
            if j != i and K_np[i, j] > 0:
                G.add_edge(cells[i], cells[j], weight=K_np[i, j])

    # ---- Community detection + layout ----
    if method == "leiden":
        node_to_idx = {n: i for i, n in enumerate(G.nodes())}
        idx_to_node = {i: n for n, i in node_to_idx.items()}
        edges = [
            (node_to_idx[u], node_to_idx[v], d["weight"])
            for u, v, d in G.edges(data=True)
        ]
        g = ig.Graph(
            edges=[(u, v) for u, v, _ in edges],
            edge_attrs={"weight": [w for _, _, w in edges]},
            directed=False,
        )
        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=res,
            seed=seed,
        )
        for idx, comm in enumerate(partition.membership):
            G.nodes[idx_to_node[idx]]["community"] = comm

        communities = [
            {n for n, d in G.nodes(data=True) if d["community"] == cid}
            for cid in sorted(set(partition.membership))
        ]

        # Community-level graph for macro layout
        G_comm = nx.Graph()
        for cid in range(len(communities)):
            G_comm.add_node(cid)
        for u, v, d in G.edges(data=True):
            cu, cv = G.nodes[u]["community"], G.nodes[v]["community"]
            if cu != cv:
                w = d.get("weight", 1.0)
                if G_comm.has_edge(cu, cv):
                    G_comm[cu][cv]["weight"] += w
                else:
                    G_comm.add_edge(cu, cv, weight=w)

        pos_comm = nx.spring_layout(
            G_comm, k=spring_k, iterations=spring_iterations_community, seed=seed
        )
        pos = {}
        for cid, nodes in enumerate(communities):
            subG = G.subgraph(nodes)
            local_pos = nx.spring_layout(
                subG, k=spring_k, iterations=spring_iterations_nodes, seed=seed
            )
            center = pos_comm[cid]
            for n, p in local_pos.items():
                pos[n] = center + compression * p

    elif method == "louvain":
        cl = CommunityLayout(
            G,
            community_compression=compression,
            layout_algorithm=nx.spring_layout,
            layout_kwargs={"k": spring_k, "iterations": spring_iterations_community},
            community_algorithm=nx.algorithms.community.louvain_communities,
            community_kwargs={"resolution": res, "seed": seed, "weight": "weight"},
        )
        pos = cl.full_positions
        communities = list(cl.communities())

    else:
        raise ValueError(f"Unknown method '{method}'. Choose 'leiden' or 'louvain'.")

    # ---- Build directed signed graph for drawing ----
    G_signed_dir = nx.DiGraph()
    for _, row in CCI_table.iterrows():
        if method == "louvain" and row["weight"] == 0:
            continue
        G_signed_dir.add_edge(row["from"], row["to"], weight=row["weight"])

    weights = np.array([d["weight"] for _, _, d in G_signed_dir.edges(data=True)])
    norm = TwoSlopeNorm(
        vmin=-abs(weights).max(), vcenter=0.0, vmax=abs(weights).max()
    )
    edge_colors = cm.coolwarm(norm(weights))
    edge_widths = 0.6 + edge_width_scale * np.abs(weights) / np.max(np.abs(weights))

    node_comm = {n: cid for cid, nodes in enumerate(communities) for n in nodes}
    n_comms = len(communities)
    cmap = get_community_cmap(n_comms)
    nodes = list(G_signed_dir.nodes())
    node_colors = [cmap(node_comm[n]) for n in nodes]
    edge_alphas = [
        edge_alpha_intra if node_comm[u] == node_comm[v] else edge_alpha_inter
        for u, v in G_signed_dir.edges()
    ]

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=figsize)

    nx.draw_networkx_nodes(
        G_signed_dir, pos, nodelist=nodes, node_color=node_colors,
        node_size=node_size, edgecolors="black", linewidths=0.7, ax=ax,
    )
    for (u, v), c, w, a in zip(G_signed_dir.edges(), edge_colors, edge_widths, edge_alphas):
        nx.draw_networkx_edges(
            G_signed_dir, pos, edgelist=[(u, v)], edge_color=[c], width=w,
            alpha=a, arrows=True, arrowstyle="-|>", arrowsize=arrow_size,
            connectionstyle="arc3,rad=0.15", ax=ax,
        )

    texts = [
        ax.text(
            x, y, node, fontsize=9, ha="center", va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.5, pad=0.2),
            zorder=10,
        )
        for node, (x, y) in pos.items()
    ]
    adjust_text(
        texts, ax=ax, expand_points=(1.2, 1.2), expand_text=(1.2, 1.2),
        force_text=0.8, force_points=0.2,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.6),
    )

    comm_handles = [
        Patch(facecolor=cmap(i), edgecolor="black", label=f"Community {i}")
        for i in range(n_comms)
    ]
    ax.legend(
        handles=comm_handles, title="Communication communities",
        bbox_to_anchor=(1.02, 1.0), loc="upper left",
    )

    sm = cm.ScalarMappable(cmap=cm.coolwarm, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.4)
    cbar.set_label("Differential LR interaction score")

    ax.set_title(
        f"Cell–cell communication community network\n"
        f"{method.capitalize()} resolution: {res}",
        fontsize=14,
    )
    ax.axis("off")
    plt.tight_layout()
    plt.show()