#!/usr/bin/env python3
import argparse
import os
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

from util import fill_phylo_ranks

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

# Format dataframe
df = pd.read_csv(os.path.join(args.i, 'IDP_perc_meta.csv'),
                 keep_default_na=False, na_values=None)
df['IDP_ratio_species'] = df['TOTAL_IDP']/df['TOTAL_PROTEIN']

fill_phylo_ranks(df)

# Scatterplot with plotly
fig = px.scatter(df, x='ASM_LENGTH',
                 y='IDP_ratio_species', color=phylo_rank,
                 labels={
                     'ASM_LENGTH': 'Assembly length (bp)',
                     'IDP_ratio_species': 'IDP ratio in species',
                     f'{phylo_rank}': phylo_rank.capitalize()
                 })

# Check output path
outpath = pathlib.Path(args.i, 'figure')
outpath.mkdir(parents=True, exist_ok=True)
fig.write_html(os.path.join(
    outpath, f'scatterplot_IDP_ratio_groupby_{phylo_rank}.html'))
