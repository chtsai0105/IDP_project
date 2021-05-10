#!/usr/bin/env python3
import argparse
import os
import pathlib
import re
from functools import partial
from math import floor, log10
from multiprocessing import Pool

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import seaborn as sns

from util import fill_phylo_ranks


def process_file(args_tuple, phylo_rank):
    filename, species, taxo = args_tuple
    tmp_df = pd.read_csv(os.path.join(args.i, 'IDP', f'{filename}.csv'))
    tmp_df['Species'] = species
    tmp_df[phylo_rank] = taxo
    return tmp_df


# Parser
avaliable_ranks = ['phylum', 'subphylum', 'class', 'subclass']
data_types = ['IDR_percentage', 'IDR_counts', 'IDR_length']
parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', metavar='<Input dir>',
                    help='Input directory (filter_{number})', required=True)
parser.add_argument('--group', '-g', metavar='<phylo group>', choices=avaliable_ranks, default='subphylum',
                    help=f'x axis catogorized by {avaliable_ranks} (default=\'subphylum\')')
parser.add_argument('--data', '-d', metavar='<datatype>', choices=data_types, default='IDR_percentage',
                    help=f'y axis measurement {data_types} (default=\'IDR_percentage\')')
parser.add_argument('--threads', '-t', metavar='', default=1, type=int,
                    help='Threads (default=1)')
args = parser.parse_args()

# Load metadata
meta_df = pd.read_csv(os.path.join(args.i, 'IDP_perc_meta.csv'),
                      keep_default_na=False, na_values=None)

fill_phylo_ranks(meta_df)
phylo_rank = args.group

# create dataframe for boxplot
df = pd.DataFrame()
args_tuples = []    # list of tuple which contains "filename", "species" and "phylo_rank"
for _, value in meta_df.iterrows():
    args_tuples.append(
        (value['FILENAME'], value['species'], value[phylo_rank]))

process_pool = Pool(processes=args.threads)
df = pd.concat(process_pool.map(
    partial(process_file, phylo_rank=phylo_rank), args_tuples))
process_pool.close()

# Grouping for event counts and average of data
grouping_df = df[[phylo_rank, args.data]].groupby(
    phylo_rank).agg(['count', 'mean'])
grouping_df.columns = ['counts', f'{args.data}_mean']
grouping_df = grouping_df.sort_values(f'{args.data}_mean')

# Merge the mean value back to df for sorting
df = df.merge(grouping_df, left_on=phylo_rank, right_index=True)
df = df.sort_values(by=f'{args.data}_mean')

# Customized colorbar
count_min = grouping_df['counts'].min()
count_max = grouping_df['counts'].max()
cmap = plt.get_cmap('Reds')
norm = mpl.colors.LogNorm(10 ** floor(log10(count_min)),
                          round(count_max * 1.1, -floor(log10(count_max))))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Customize color palette by IDR numbers
customized_palette = grouping_df['counts'].to_dict()
customized_palette = dict(
    map(lambda x: (x[0], cmap(norm(x[1]))), customized_palette.items()))

# Plot
plt.figure(figsize=(15, 8))
ax = sns.boxplot(data=df, x=phylo_rank, y=args.data,
                 showfliers=False, palette=customized_palette)
plt.xlabel(phylo_rank.capitalize())
plt.ylabel(re.sub('_', ' ', args.data))

# Reformat xticklabels with more info if x axis items fewer than 20
if len(grouping_df) <= 20:
    text_lis = []
    for t in ax.get_xticklabels():
        taxa = t.get_text()
        sp_num = len(df.loc[df[phylo_rank] == taxa, 'Species'].unique())
        idr_num = grouping_df['counts'][taxa]
        new_text = f'{taxa}\n({idr_num:,} idr among {sp_num:,} species)'
        text_lis.append(new_text)
    ax.set_xticklabels(text_lis)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")

# Add colorbar to represent number of IDRs
cbar = plt.colorbar(sm)
cbar.set_label('# of IDRs', rotation=270, labelpad=15)

plt.tight_layout()

# Output figure
outpath = pathlib.Path(args.i, 'figure')
outpath.mkdir(parents=True, exist_ok=True)
plt.savefig(os.path.join(outpath, f'boxplot_{phylo_rank}-{args.data}.png'))
