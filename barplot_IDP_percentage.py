#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib as mpl
import seaborn as sns
import re


def scale_range(input, min=0, max=1):
    norm = input - (np.min(input))
    norm /= np.max(norm) / (max - min)
    norm += min
    return norm


df = pd.read_csv('./lib/IDP_percentage.csv')
df['IDP_percentage'] = df['Total_IDP']/df['Total_protein']
IDP_mean = df.groupby('Species').mean()
IDP_mean.rename(columns={"IDP_percentage": "IDP_mean"}, inplace=True)
df = df.merge(IDP_mean['IDP_mean'], left_on='Species',
              right_index=True).sort_values(by='IDP_mean')
start = 0
end = 0
lis = []
start_idx = []
end_idx = []
for species in df['Species']:
    if species not in lis:
        lis.append(species)
    end += 1
    if len(lis) == 51:
        start_idx.append(start)
        end_idx.append(end - 1)
        start = end
        lis = []
    elif end == len(df['Species']):
        start_idx.append(start)
        end_idx.append(end)

col_pal = sns.color_palette("rocket_r", n_colors=256)
color_idx = scale_range(df['IDP_mean'].values, max=255).round().astype(int)
cmap = plt.get_cmap('rocket_r')
norm = mpl.colors.Normalize(
    vmin=df['IDP_mean'].min()*100, vmax=df['IDP_mean'].max()*100)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

for idx, (start, end) in enumerate(zip(start_idx, end_idx)):
    plt.figure(figsize=(8, 10))
    cbar = plt.colorbar(sm, format='%.0f%%')
    ax = sns.barplot(y='Species', x='IDP_percentage',
                     data=df[start: end], palette=col_pal[color_idx[start]: color_idx[end]])
    plt.xlim(0, 1)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(1))
    sns.despine()
    plt.savefig(f'./figure/IDP_protein_ratio_{idx}.png', bbox_inches='tight')
