def fill_phylo_ranks(df):
    desired_ranks = ['phylum', 'subphylum', 'class',
                     'subclass', 'order', 'family', 'genus', 'species']
    # Fill empty ranks with higher class names
    for col in df.loc[:, desired_ranks]:
        na_idx = df[col] == ""
        if na_idx.any() == True:
            df.loc[na_idx, col] = df[previous_col][na_idx]
        previous_col = col
