import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import *
from adjustText import adjust_text

def plot_cci(graph, colors, plt_name, coords, pg, emax=None, leg=False, low=25, high=75, ignore_alpha=False, log=False, efactor=8, vfactor=12, vnames=True,figsize=None):
    
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
    node_sizes = [size*500 for size in pg]

    # Plot the graph
    if figsize is None:
        figsize = (7,7)
    fig, ax = plt.subplots(figsize=figsize)
    
    nx.draw(graph, pos=coords_scale, edge_color=edge_colors, node_color=node_colors, node_size=node_sizes,
            width=[d['width'] for _, _, d in graph.edges(data=True)],
            arrows=True, arrowsize=15, arrowstyle='-|>',
            connectionstyle='arc3,rad=0.3', with_labels=vnames)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)

    # Node Pagerank legend
    if pg is not None:
        min_pg, max_pg = min(pg), max(pg)
        legend1 = ax.legend(loc='lower left', title="Pagerank",
                handles=[plt.Line2D([], [], linestyle='', marker='o', markersize=v / vfactor, markerfacecolor='black', markeredgecolor='none') for v in [min_pg, (min_pg + max_pg) / 2, max_pg]],
                labels=[round(min_pg, 2), round((min_pg + max_pg) / 2, 2), round(max_pg, 2)],  bbox_to_anchor=(0.74, 0.5))

    # Thickness legend
    non_zero_inter_edges = [d['inter'] for _, _, d in graph.edges(data=True) if d.get('inter', 0) != 0]
    if non_zero_inter_edges:
        e_wid_sp = [round(min(non_zero_inter_edges), 2), round(min(non_zero_inter_edges) + (emax / 2), 2), round(emax, 2)]
        legend2 = ax.legend(e_wid_sp, title='Percentage of \nthe interactions', title_fontsize='small', loc='upper left', bbox_to_anchor=(0.74, 0.5))

    ax.add_artist(legend1)
    ax.add_artist(legend2)
    
    ax.set_title(plt_name)
    # Show the plot
    plt.tight_layout()
    plt.show()

def plot_pca_LR_comparative(lrobj_tblPCA, pca_table, dims=(1, 2), ret=False, ggi=True, include_tf=False, gene_types="all"):
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