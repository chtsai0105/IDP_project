#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns
import plotly.express as px


# Merge metadata and IDP table
meta_df = pd.read_csv('./lib/ensembl_metadata.csv', keep_default_na=False, na_values=None)
IDP_df = pd.read_csv('./lib/IDP_percentage.csv')

merged_df = pd.merge(meta_df, IDP_df, on='FILENAME')
merged_df.to_csv('./lib/meta_IDP_perc.csv', index=False)

# Format dataframe
merged_df = pd.read_csv('./lib/meta_IDP_perc.csv', keep_default_na=False, na_values=None)
merged_df['IDP_perc'] = merged_df['TOTAL_IDP'] / merged_df['TOTAL_PROTEIN']
phylo_df = merged_df[['NAME', 'ASM_LENGTH', 'IDP_perc', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus']]

# Fill empty columns with higher level taxonomic ranks
for col in phylo_df.iloc[:, 3:]:
    na_idx = phylo_df[col] == ""
    if na_idx.any() == True:
        phylo_df.loc[na_idx, col] = phylo_df[previous_col][na_idx]
    previous_col = col

# Scatterplot with plotly
fig = px.scatter(phylo_df, x="ASM_LENGTH", y="IDP_perc", color="subphylum")
fig.write_html("./figure/")
