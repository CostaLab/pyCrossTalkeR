from Python import plot as plot
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load data from pickle file
with open("output/LR_data.pkl", "rb") as f:
    data = pickle.load(f)

# Cell-Cell-Interaction Analysis
plot.plot_cci(graph=data["graphs"]["EXP"],
         colors=data["colors"],
         plt_name='Disease',
         coords=data["coords"],
         emax= None,
         leg= False,
         low= 0,
         high= 0, 
         ignore_alpha= False,
         log= False,
         efactor= 6,
         vfactor= 0.03, 
         pg= data["rankings"]["EXP"]["Pagerank"]
         )

plot.plot_cci(graph=data["graphs"]["CTR"],
         colors=data["colors"],
         plt_name='Control',
         coords=data["coords"],
         emax= None,
         leg= False,
         low= 0,
         high= 0, 
         ignore_alpha= False,
         log= False,
         efactor= 6,
         vfactor= 0.03, 
         pg= data["rankings"]["CTR"]["Pagerank"]
         )

# Analysis on the Gene-Cell Interaction Level

rankings_table_ctr = data['rankings']['CTR_ggi'].sort_values(by='Pagerank', ascending=False).head(10)
rankings_table_ctr['signal'] = ['negative' if x < 0 else 'positive' for x in rankings_table_ctr['Pagerank']]
rankings_table_exp = data['rankings']['EXP_ggi'].sort_values(by='Pagerank', ascending=False).head(10)
rankings_table_exp['signal'] = ['negative' if x < 0 else 'positive' for x in rankings_table_exp['Pagerank']]
custom_palette = {'positive': '#ff3e3e', 'negative': '#00FFFF'}  # Orange and Blue
fig, axs = plt.subplots(1, 2, figsize=(15, 8))

# Plot for CTR_ggi
sns.barplot(ax=axs[0], x='Pagerank', y='nodes', data=rankings_table_ctr, hue='signal', dodge=False, palette=custom_palette)
axs[0].set_title("Top Listener in Control Condition")
axs[0].set_xlabel('Pagerank')
axs[0].set_ylabel('Nodes')

# Set x-axis tick intervals for CTR_ggi
ctr_max = rankings_table_ctr['Pagerank'].max()
ctr_min = rankings_table_ctr['Pagerank'].min()
ctr_ticks = np.linspace(ctr_min, ctr_max, num=3)  # Adjust 'num' for more/less intervals
axs[0].set_xticks(ctr_ticks)
axs[0].set_xticklabels([f'{tick:.2f}' for tick in ctr_ticks])

# Plot for EXP_ggi
sns.barplot(ax=axs[1], x='Pagerank', y='nodes', data=rankings_table_exp, hue='signal', dodge=False, palette=custom_palette)
axs[1].set_title("Top Listener in Disease Condition")
axs[1].set_xlabel('Pagerank')
axs[1].set_ylabel('')

# Set x-axis tick intervals for EXP_ggi
exp_max = rankings_table_exp['Pagerank'].max()
exp_min = rankings_table_exp['Pagerank'].min()
exp_ticks = np.linspace(exp_min, exp_max, num=3)  # Adjust 'num' for more/less intervals
axs[1].set_xticks(exp_ticks)
axs[1].set_xticklabels([f'{tick:.2f}' for tick in exp_ticks])

axs[0].grid(True, linestyle='--', linewidth=0.5)
axs[1].grid(True, linestyle='--', linewidth=0.5)
axs[0].set_axisbelow(True)
axs[1].set_axisbelow(True)
handles, labels = axs[0].get_legend_handles_labels()
plt.tight_layout()
plt.show()