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
species_dict = dict()
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
    file_name = '.'.join((url_name, asm_def))
    # URL for fasta file
    dbname = re.sub(rf'_core_{release}_\d+_\d+', "",
                    dat['databases'][0]['dbname'])
    subfolder = dat['organism']['name']
    # Handle collection folders
    if dbname.startswith('fungi_'):
        subfolder = '/'.join((dbname, subfolder))
    # Final url (Unfinished, currently point to the folder but not fasta file yet)
    dna_url = '/'.join((url_base, 'fasta', subfolder, 'dna'))
    pep_url = '/'.join((url_base, 'fasta', subfolder, 'pep'))
    # Output to dict
    species_dict[idx] = [asm_acc, display_name, asm_name, taxo_id,
                         asm_length, total_prot, dna_url, pep_url, file_name]

df = pd.DataFrame.from_dict(species_dict, orient='index', columns=[
                            'ACCESSION', 'SPECIES', 'ASM_NAME', 'NCBI_TAXID', 'ASM_LENGTH', 'TOTAL_PROTEIN', 'DNA_URL', 'PEP_URL', 'FILENAME'])


ncbi = NCBITaxa()
taxids = df['NCBI_TAXID']
desired_ranks = ['phylum', 'subphylum', 'class',
                 'subclass', 'order', 'family', 'genus']

results = list()
for taxid in taxids:
    results.append(get_desired_ranks(taxid, desired_ranks))
taxo_df = pd.DataFrame(results, columns=desired_ranks)

df = pd.concat([df, taxo_df], axis=1)
df.to_csv('./lib/ensembl_metadata.csv', index=False)
