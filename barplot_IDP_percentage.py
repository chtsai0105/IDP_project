#!/usr/bin/env python3
import argparse
import re
import os
import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import seaborn as sns

from util import fill_phylo_ranks


def scale_range(input, min=0, max=1):
    norm = input - (np.min(input))
    norm /= np.max(norm) / (max - min)
    norm += min
    return norm


# Parser
avaliable_ranks = ['phylum', 'subphylum', 'class',
                   'subclass', 'order', 'family', 'genus', 'species']
parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', metavar='<Input dir>',
                    help='Input directory (filter_{number})', required=True)
parser.add_argument('--group', '-g', metavar='<phylo group>', choices=avaliable_ranks, default='subphylum',
                    help=f'Used for catogorized y axis {avaliable_ranks} (default=\'subphylum\')')
args = parser.parse_args()
phylo_rank = args.group

# Generate dataframe for plot
df = pd.read_csv(os.path.join(args.i, 'IDP_perc_meta.csv'),
                 keep_default_na=False, na_values=None)
df['IDP_ratio_species'] = df['TOTAL_IDP']/df['TOTAL_PROTEIN']

fill_phylo_ranks(df)

IDP_mean = df.groupby(phylo_rank).mean()
IDP_mean.rename(columns={"IDP_ratio_species": "IDP_ratio_mean"}, inplace=True)
df = df.merge(IDP_mean['IDP_ratio_mean'], left_on=phylo_rank,
              right_index=True).sort_values(by='IDP_ratio_mean')

# Set indeces for 50 entries per figure
start = 0
end = 0
lis = []
start_idx = []
end_idx = []
for species in df[phylo_rank]:
    if species not in lis:
        lis.append(species)
    end += 1
    if len(lis) == 51:
        start_idx.append(start)
        end_idx.append(end - 1)
        start = end
        lis = []
    elif end == len(df[phylo_rank]):
        start_idx.append(start)
        end_idx.append(end)

# Set colormap
cmap = plt.get_cmap('Reds')
norm = mpl.colors.Normalize(vmin=0, vmax=df['IDP_ratio_mean'].max() * 100)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Set customized color palette
customized_palette = {}
for sp, value in IDP_mean.iterrows():
    customized_palette[sp] = cmap(norm(value['IDP_ratio_mean'] * 100))

# Check output path
outpath = pathlib.Path(args.i, 'figure')
outpath.mkdir(parents=True, exist_ok=True)

# Plot
for idx, (start, end) in enumerate(zip(start_idx, end_idx)):
    plt.figure(figsize=(8, 10))
    cbar = plt.colorbar(sm, format='%.0f%%')
    ax = sns.barplot(y=phylo_rank, x='IDP_ratio_species',
                     data=df[start: end], palette=customized_palette)
    plt.xlim(0, 1)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(1))
    plt.xlabel(f'IDP ratio in {phylo_rank}')
    plt.ylabel(phylo_rank.capitalize())
    sns.despine()
    if len(start_idx) == 1:
        idx = ""
    else:
        idx = f'_{idx}'
    plt.savefig(os.path.join(
        outpath, f'barplot_IDP_ratio_in_{phylo_rank}{idx}.png'), bbox_inches='tight')
