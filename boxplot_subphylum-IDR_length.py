#!/usr/bin/env python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

meta_df = pd.read_csv('./lib/ensembl_metadata.csv',
                      keep_default_na=False, na_values=None)
for col in meta_df.iloc[:, 9:]:
    na_idx = meta_df[col] == ""
    if na_idx.any() == True:
        meta_df.loc[na_idx, col] = meta_df[previous_col][na_idx]
    previous_col = col

df = pd.DataFrame()
for idx, filename in enumerate(meta_df['FILENAME']):
    tmp_df = pd.read_csv(f'./output/{filename}_IDR.csv')
    tmp_df['IDR_length'] = tmp_df['IDR_end'] - tmp_df['IDR_start']
    tmp_df = tmp_df.groupby('IDP').sum().reset_index()[
        ['IDR_length', 'IDR_percentage']]
    tmp_df['Species'] = meta_df.loc[idx, 'SPECIES']
    tmp_df['subphylum'] = meta_df.loc[idx, 'subphylum']
    df = pd.concat([df, tmp_df])

# Customized colorbar
cmap = plt.get_cmap('Reds')
norm = mpl.colors.LogNorm(1000, 3000000)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Customize color palette by IDR numbers
customized_palette = {}
for sp in df['subphylum'].unique():
    idr_num = len(df[df['subphylum'] == sp])
    customized_palette[sp] = cmap(norm(idr_num))

plt.figure(figsize=(15, 8))
ax = sns.boxplot(data=df, x='subphylum', y='IDR_length',
                 showfliers=False, palette=customized_palette)

# Reformat xticklabels
text_lis = []
for t in ax.get_xticklabels():
    subphylum = t.get_text()
    sp_num = len(df.loc[df['subphylum'] == subphylum, 'Species'].unique())
    idr_num = len(df[df['subphylum'] == subphylum])
    new_text = f'{subphylum}\n({idr_num:,} IDPs among {sp_num:,} species)'
    text_lis.append(new_text)
ax.set_xticklabels(text_lis, rotation=30, ha="right")

# Add colorbar to represent number of IDRs
cbar = plt.colorbar(sm)
cbar.set_label('# of IDRs', rotation=270, labelpad=15)

plt.tight_layout()
plt.savefig('./figure/boxplot_subphylum-IDR_length.png')
