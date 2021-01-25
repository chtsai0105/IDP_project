#!/usr/bin/env python3
import gzip

import iupred2a_lib
import numpy as np
import pandas as pd
from Bio import SeqIO


def idr_detector(arr, thresh=0.5, window_size=30):
    m = arr >= thresh
    me = np.r_[False, m, False]
    idx = np.flatnonzero(me[:-1] != me[1:])
    regions = idx.reshape((int(idx.shape[0]/2), 2))
    lens = regions[:, 1] - regions[:, 0]
    return regions[(lens >= window_size)]

# Metadata, would be constructed with dataframe format in the future
species_lis = ['Neurospora_crassa', 'Aspergillus_fumigatus']
fasta_lis = ['pep/Neurospora_crassa.NC12.pep.all.fa.gz', 'pep/Aspergillus_fumigatus.ASM265v1.pep.all.fa.gz']

df = pd.DataFrame(columns=['IDP', 'IDR_percentage', 'IDR', 'Species'])
for species, fasta in zip(species_lis, fasta_lis):
    sub_df = pd.DataFrame(columns=['IDP', 'IDR_percentage', 'IDR'])
    with gzip.open(fasta, 'rt') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            iupred_scores = np.asarray(iupred2a_lib.iupred(record.seq, "long"))
            IDR_in_entry = idr_detector(iupred_scores)
            if IDR_in_entry.size != 0:
                IDR_lens = IDR_in_entry[:, 1] - IDR_in_entry[:, 0]
                IDR_percentage = IDR_lens.sum() / len(record.seq)
                lis = IDR_in_entry.tolist()
                string = ';'.join(str(e) for e in lis)
                string = string.replace(',', ':')
                tmp_df = pd.Series([record.id, IDR_percentage, string], index=['IDP', 'IDR_percentage', 'IDR'])
                sub_df = pd.concat([sub_df, tmp_df], ignore_index=True)
        sub_df['Species'] = species
    df = pd.concat([df, sub_df], ignore_index=True)

# Reorder the dataframe
df = df[['Species', 'IDP', 'IDR_percentage', 'IDR']]

# Save to csv file
df.to_csv('result_IDR_in_column.csv', index=False)
