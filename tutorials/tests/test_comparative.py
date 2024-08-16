
from Python import plot as plot
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load data from pickle file
with open("output/LR_data.pkl", "rb") as f:
    data = pickle.load(f)


# Cell-Cell-Interaction analysis results
plot.plot_cci(graph=data["graphs"]["EXP_x_CTR"],
         colors=data["colors"],
         plt_name='Disease vs Control',
         coords=data["coords"],
         emax= None,
         leg= False,
         low= 0,
         high= 0, 
         ignore_alpha= False,
         log= False,
         efactor= 6,
         vfactor= 0.03, 
         pg= data["rankings"]["EXP_x_CTR"]["Pagerank"]
         )


custom_palette = {'positive': '#00FFFF', 'negative':'#ff3e3e'}  # Orange and Blue

# Iterate through the rankings and plot
for key in data['rankings']:
    if '_x_' in key and 'ggi' not in key:
        rankings_table = data['rankings'][key]
        rankings_table = rankings_table.sort_values(by='Pagerank')
        rankings_table['signal'] = ['negative' if x < 0 else 'positive' for x in rankings_table['Pagerank']]

        # Plot
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Pagerank', y='nodes', data=rankings_table, hue='signal', dodge=False, palette=custom_palette)
        plt.title(f"Ranking for {key}")
        plt.xlabel('Pagerank')
        plt.ylabel('Nodes')

        # Set x-axis tick intervals
        max_val = rankings_table['Pagerank'].max()
        min_val = rankings_table['Pagerank'].min()
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
        


for key in data['rankings']:
    if '_x_' in key and 'ggi' not in key:
        rankings_table = data['rankings'][key]
        rankings_table = rankings_table.sort_values(by='Influencer')
        rankings_table['signal'] = ['negative' if x < 0 else 'positive' for x in rankings_table['Influencer']]

        plt.figure(figsize=(10, 6))
        sns.barplot(x='Influencer', y='nodes', data=rankings_table, hue= 'signal', dodge= False, palette = custom_palette)
        

        plt.title(key)
        plt.xlabel('Influencer')
        plt.ylabel('Nodes')

        # Set x-axis tick intervals
        max_val = rankings_table['Influencer'].max()
        min_val = rankings_table['Influencer'].min()
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


#Cell-Gene-Interaction Analysis Results

# Iterate through the rankings and plot
for key in data['rankings']:
    if '_x_' in key and 'ggi' in key:
        rankings_table = data['rankings'][key]
        rankings_table = rankings_table.loc[rankings_table['Pagerank'].abs().nlargest(20).index]
        rankings_table = rankings_table.sort_values(by='Pagerank')
        rankings_table['signal'] = ['negative' if x < 0 else 'positive' for x in rankings_table['Pagerank']]

        # Plot
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Pagerank', y='nodes', data=rankings_table, hue='signal', dodge=False, palette=custom_palette)
        plt.title(f"Ranking for {key}")
        plt.xlabel('Pagerank')
        plt.ylabel('Nodes')

        # Set x-axis tick intervals
        max_val = rankings_table['Pagerank'].max()
        min_val = rankings_table['Pagerank'].min()
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

plot.plot_pca_LR_comparative(
    lrobj_tblPCA = data,
    pca_table = "EXP_x_CTR_ggi",
    dims = (1, 2),
    ret = True,
    ggi = True
)
