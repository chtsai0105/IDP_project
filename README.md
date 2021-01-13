# IDP_project
Intrinsic disordered proteins (IDP) in fungi project in Stajich lab

S288c ORF translations was downloaded from [SGD](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/)

## Requirement
Biopython(>=1.76)
python(>=3.6.7)

## Pipeline
All the analysis is done on HPCC. First, load the IUPred module

```shell
module load iupred
```

Run the python script
```shell
bash run_iupred_genome-wide.py
```

A region that have at least 30 consecutive residues with iupred score over 0.5 would be determined as putative disordered region (IDR). 
Any protein that have at least 1 IDR would be classified as IDP. 
For all IDPs, the name of the protein, IDR percentage among all residues as well the start and end residue of each IDR would be reported and stored in a pickle file.
Finally, the script would conclude how many proteins in these proteome would possibly be IDP and the average disordered region length of all IDP.
