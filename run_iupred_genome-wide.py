#!/usr/bin/env python3
import gzip
import pickle

import iupred2a_lib
import numpy as np
from Bio import SeqIO

import numpy as np
import pandas as pd
from util.py import connect_to_db

def idr_detector(arr, thresh=0.5, window_size=30):
    m = arr >= thresh
    me = np.r_[False, m, False]
    idx = np.flatnonzero(me[:-1] != me[1:])
    regions = idx.reshape((int(idx.shape[0]/2), 2))
    lens = regions[:, 1] - regions[:, 0]
    return regions[(lens >= window_size)]

species = 'Neurospora_crassa'
fasta = 'pep/Neurospora_crassa.NC12.pep.all.fa.gz'

df = pd.DataFrame(columns=['IDR_start', 'IDR_end', 'IDP'])
with gzip.open(fasta, 'rt') as fh:
    for record in SeqIO.parse(fh, "fasta"):
        iupred_scores = np.asarray(iupred2a_lib.iupred(record.seq, "long"))
        IDR_in_entry = idr_detector(iupred_scores)
        if IDR_in_entry.size != 0:
            IDR_lens = IDR_in_entry[:, 1] - IDR_in_entry[:, 0]
            # IDR_percentage.append(IDR_lens.sum() / len(record.seq))
            tmp_df = pd.DataFrame(IDR_in_entry, columns=['IDR_start', 'IDR_end'])
            tmp_df['IDP'] = record.id
            df = df.append(tmp_df, ignore_index=True)

# Reorder the dataframe
df = df[['IDP', 'IDR_start', 'IDR_end']]

# Save to database
engine = connect_to_db()
df.to_sql(species, con=engine, if_exists='replace', index_label='ID')
