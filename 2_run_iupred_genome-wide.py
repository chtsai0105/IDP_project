#!/usr/bin/env python3
import argparse
import gzip
import os
import pathlib
import re

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
parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', metavar='<Input.gz>',
                    help='Peptide fasta.gz', required=True)
parser.add_argument('--o', metavar='<Output dir>',
                    help='Directory of output csv', required=True)
args = parser.parse_args()

# Remove prefix and postfix in input path
fasta = args.i
species = os.path.basename(fasta)
species = re.sub('.pep.all.fa.gz$', '', species)


idr_df = pd.DataFrame()
with gzip.open(fasta, 'rt') as fh:
    for record in SeqIO.parse(fh, "fasta"):
        iupred_scores = np.asarray(iupred2a_lib.iupred(record.seq, "long"))
        IDR_in_entry = idr_detector(iupred_scores)
        if IDR_in_entry.size != 0:
            IDR_lens = IDR_in_entry[:, 1] - IDR_in_entry[:, 0]
            IDR_percentage = IDR_lens / len(record.seq)
            tmp_df = pd.DataFrame(IDR_in_entry, columns=[
                                  'IDR_start', 'IDR_end'])
            tmp_df['protein_length'] = len(record.seq)
            tmp_df['IDR_length'] = IDR_lens
            tmp_df['IDR_percentage'] = IDR_percentage
            tmp_df['IDP'] = record.id
            idr_df = pd.concat([idr_df, tmp_df], ignore_index=True)

# Reorder the dataframe
idr_df = idr_df[['IDP', 'protein_length', 'IDR_start',
                 'IDR_end', 'IDR_length', 'IDR_percentage']]

# Save to csv file
pathlib.Path(args.o).mkdir(parents=True, exist_ok=True)
idr_df.to_csv(os.path.join(args.o, f'{species}.csv'), index=False)
