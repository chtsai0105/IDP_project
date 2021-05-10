#!/usr/bin/env python

import argparse
import os
import pathlib
import re
from functools import partial
from multiprocessing import Pool

import pandas as pd
import psutil


def get_files(directory, pattern):
    for path in pathlib.Path(directory).rglob(pattern):
        yield path.absolute()


def process_file(filename, out_path, window=30):
    df = pd.read_csv(filename)
    df = df[df['IDR_length'] >= window]

    grouping_df = df.groupby('IDP').agg(['count', 'sum', 'mean'])
    grouping_df = grouping_df.loc[:, [('protein_length', 'mean'), (
        'IDR_length', 'count'), ('IDR_length', 'sum'), ('IDR_percentage', 'sum')]]
    grouping_df.columns = ['protein_length',
                           'IDR_counts', 'IDR_length', 'IDR_percentage']
    grouping_df = grouping_df.reset_index()

    species = os.path.basename(filename)
    species = re.sub('.csv$', '', species)
    total_IDP = len(grouping_df)

    # Output IDP.csv
    grouping_df.to_csv(os.path.join(out_path, f'{species}.csv'), index=False)
    return [species, total_IDP]


# Parser
parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', metavar='<Input dir>',
                    help='Directory of input (IDR) csv', required=True)
parser.add_argument('--o', metavar='<Output dir>',
                    help='Output directory', required=True)
parser.add_argument('-w', '--window', metavar='', default=30, type=int,
                    help='Minimum number of residues to define IDR (default=30)')
parser.add_argument('-t', '--threads', metavar='', default=1, type=int,
                    help='Threads (default=1)')
args = parser.parse_args()

# Check whether output path exist
IDPcsv_outpath = pathlib.Path(args.o, 'IDP')
IDPcsv_outpath.mkdir(parents=True, exist_ok=True)

# Get files based on pattern
files = []
for file in get_files(directory=args.i, pattern='*.csv'):
    files.append(file)

# Multiprocess
process_pool = Pool(processes=args.threads)
data = process_pool.map(
    partial(process_file, out_path=IDPcsv_outpath, window=args.window), files)
process_pool.close()

# IDP table
IDP_perc_df = pd.DataFrame(data, columns=['FILENAME', 'TOTAL_IDP'])

# Merge metadata and IDP table
meta_df = pd.read_csv('./lib/ensembl_metadata.csv',
                      keep_default_na=False, na_values=None)

merged_df = pd.merge(meta_df, IDP_perc_df, on='FILENAME')
merged_df.to_csv(os.path.join(
    args.o, f'IDP_perc_meta.csv'), index=False)
