#!/usr/bin/env python3
from Bio import SeqIO
import iupred2a_lib
import pickle
import numpy as np

def idr_detector(arr, thresh=0.5, window_size=30):
    m = arr >= thresh
    me = np.r_[False, m, False]
    idx = np.flatnonzero(me[:-1] != me[1:])
    regions = idx.reshape((int(idx.shape[0]/2), 2))
    lens = regions[:, 1] - regions[:, 0]
    return regions[(lens >= window_size)]

inputfile = "orf_trans.fasta"
IDP = []
IDR_percentage = []
IDR = []

for record in SeqIO.parse(inputfile, "fasta"):
    iupred_scores = np.asarray(iupred2a_lib.iupred(record.seq, "long"))
    IDR_in_entry = idr_detector(iupred_scores)
    if IDR_in_entry.size != 0:
        IDP.append(record.id)
        IDR_lens = IDR_in_entry[:, 1] - IDR_in_entry[:, 0]
        IDR_percentage.append(IDR_lens.sum() / len(record.seq))
        IDR.append(IDR_in_entry)

# Save output to pickle file
with open('IDP_S288c.pickle', 'wb') as fh:
    pickle.dump((IDP, IDR_percentage, IDR), fh)

parser = list(SeqIO.parse(inputfile, "fasta"))
lens_of_total_protein = len(parser)
print("Disordered protein : {} out of {}".format(len(IDR_percentage, lens_of_total_protein)))
print("Average disordered regions length of IDP: {}".format(np.average(np.asarray(IDR_percentage))))
