import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
import plotly.colors as pc
import plotly.graph_objects as go
from plotnine import *
from adjustText import adjust_text
from gprofiler import GProfiler
from sankeyflow import Sankey


def plot_cci(graph, colors, plt_name, coords, pg, emax=None, leg=False, low=25, high=75, ignore_alpha=False, log=False, efactor=8, vfactor=12, vnames=True, figsize=None, scale_factor=2, node_size=2, font_size=10):
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
        not include transparency on the plot
    log :
        logscale the interactions
    efactor :
        edge scale factor
    vfactor :
        edge scale factor
    vnames :
        remove vertex labels
    pg :
        pagerank values

    Returns
    -------
    Python default plot
    
    """

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
                labels=[round(min_pg, 2), round((min_pg + max_pg) / 2, 2), round(max_pg, 2)],  bbox_to_anchor=(0.8, 0))

    # Thickness legend
    non_zero_inter_edges = [d['inter'] for _, _, d in graph.edges(data=True) if d.get('inter', 0) != 0]
    if non_zero_inter_edges:
        e_wid_sp = [round(min(non_zero_inter_edges), 2), round(min(non_zero_inter_edges) + (emax / 2), 2), round(emax, 2)]
        legend2 = ax.legend(e_wid_sp, title='Percentage of \nthe interactions', title_fontsize='small', loc='upper left', bbox_to_anchor=(0.8, 0.4))

    ax.add_artist(legend1)
    ax.add_artist(legend2)
    
    ax.set_title(plt_name)
    # Show the plot
    plt.tight_layout()
    plt.show()


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
    

def plot_bar_rankings(data, table_name, ranking, type = None, filter_sign = None, mode = "cci", top_num = 10):
    """
    This function generates the barplot for a given network ranking on the CGI level. Further, the genes can be filtered by selected gene types to filter the plot.

    Parameters
    ----------
    data_object :
        LRobject with all data

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
        rankings_table = data['rankings'][table_name]

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


def plot_clustermap(data, title, annot=True):
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

    Returns
    -------
    Python Cluster map

    """
    pivot_table = data.groupby(["source", "target"])["LRScore"].sum().unstack().fillna(0)
    xlabel, ylabel = "Target Tissue", "Source Tissue"
    
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