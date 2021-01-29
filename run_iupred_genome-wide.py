#!/usr/bin/env python3
import gzip
import argparse

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


# Parser
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--i', metavar='<Input.gz>', help='Peptide fasta.gz')
parser.add_argument('--o', metavar='<Output.csv>', help='File name of output csv')
args = parser.parse_args()

# Remove prefix and postfix in input path
fasta = args.i
species = fasta.strip('pep/').strip('.pep.all.fa.gz')


df = pd.DataFrame(columns=['IDR_start', 'IDR_end', 'IDR_percentage', 'IDP'])
with gzip.open(fasta, 'rt') as fh:
    for record in SeqIO.parse(fh, "fasta"):
        iupred_scores = np.asarray(iupred2a_lib.iupred(record.seq, "long"))
        IDR_in_entry = idr_detector(iupred_scores)
        if IDR_in_entry.size != 0:
            IDR_lens = IDR_in_entry[:, 1] - IDR_in_entry[:, 0]
            IDR_percentage = IDR_lens / len(record.seq)
            tmp_df = pd.DataFrame(IDR_in_entry, columns=['IDR_start', 'IDR_end'])
            tmp_df['IDR_percentage'] = IDR_percentage
            tmp_df['IDP'] = record.id
            df = pd.concat([df, tmp_df], ignore_index=True)
df['Species'] = species

# Reorder the dataframe
df = df[['Species', 'IDP', 'IDR_start', 'IDR_end', 'IDR_percentage']]

# Save to csv file
df.to_csv(f'{args.o}.csv', index=False)
