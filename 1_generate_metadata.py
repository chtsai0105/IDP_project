#!/usr/bin/env python3

import json
import re

import numpy as np
import pandas as pd
from ete3 import NCBITaxa


def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid)
                         for (taxid, rank) in lineage2ranks.items())
    result = list()
    for rank in desired_ranks:
        taxid = ranks2lineage.get(rank, np.nan)
        if taxid is not np.nan:
            result.append(list(ncbi.get_taxid_translator([taxid]).values())[0])
        else:
            result.append(taxid)
    return(result)


with open('./lib/species_metadata_EnsemblFungi.json') as jsonin:
    data = json.load(jsonin)

# Setup host and release
host = 'http://ftp.ebi.ac.uk/ensemblgenomes'
release = '49'
basefolder = f'pub/release-{release}/fungi'
url_base = '/'.join((host, basefolder))

# Retrieve data from ensembl metadata
meta_dict = dict()
download_dict = dict()
for idx, dat in enumerate(data):
    # Assembly info
    display_name = dat['organism']['display_name']
    url_name = dat['organism']['url_name']
    taxo_id = dat['organism']['species_taxonomy_id']
    asm_acc = dat['assembly']['assembly_accession']
    asm_length = dat['assembly']['base_count']
    asm_name = dat['assembly']['assembly_name']
    asm_def = dat['assembly']['assembly_default']
    total_prot = dat['annotations']['nProteinCoding']
    # URL for DNA, pep fasta and gff file
    basename = '.'.join((url_name, asm_def))
    dna_urlname = '.'.join((basename, 'dna', 'toplevel', 'fa', 'gz'))
    pep_urlname = '.'.join((basename, 'pep', 'all', 'fa', 'gz'))
    gff_urlname = '.'.join((basename, release, 'gff3', 'gz'))
    # Handle subfolder name
    dbname = re.sub(rf'_core_{release}_\d+_\d+', "",
                    dat['databases'][0]['dbname'])
    subfolder = dat['organism']['name']
    # Handle collection folders
    if dbname.startswith('fungi_'):
        subfolder = '/'.join((dbname, subfolder))
    # Final url
    dna_url = '/'.join((url_base, 'fasta', subfolder, 'dna', dna_urlname))
    pep_url = '/'.join((url_base, 'fasta', subfolder, 'pep', pep_urlname))
    gff_url = '/'.join((url_base, 'gff3', subfolder, gff_urlname))
    # Output to dict
    meta_dict[idx] = [asm_acc, display_name, basename, asm_name, asm_length, total_prot, taxo_id]
    download_dict[idx] = [asm_acc, display_name, dna_url, pep_url, gff_url]

# Convert dictionary to dataframe
meta_df = pd.DataFrame.from_dict(meta_dict, orient='index', columns=[
    'ACCESSION', 'DISPLAY_NAME', 'FILENAME', 'ASM_NAME', 'ASM_LENGTH', 'TOTAL_PROTEIN', 'NCBI_TAXID'])
download_df = pd.DataFrame.from_dict(download_dict, orient='index', columns=[
                                     'ACCESSION', 'DISPLAY_NAME', 'DNA_URL', 'PEP_URL', 'GFF_URL'])
# Save the download link table
download_df.to_csv('./lib/ensembl_download.csv', index=False)

# Using ete3 package to handle taxa info
ncbi = NCBITaxa()
taxids = meta_df['NCBI_TAXID']
desired_ranks = ['phylum', 'subphylum', 'class',
                 'subclass', 'order', 'family', 'genus', 'species']

results = list()
for taxid in taxids:
    results.append(get_desired_ranks(taxid, desired_ranks))
taxo_df = pd.DataFrame(results, columns=desired_ranks)

# Join taxo_df as additional columns
meta_df = pd.concat([meta_df, taxo_df], axis=1)
# Save the metadata table
meta_df.to_csv('./lib/ensembl_metadata.csv', index=False)
