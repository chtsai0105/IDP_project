#!/bin/bash

echo "Species,Total_protein,Total_IDP" > stat/IDP_percentage.csv

for species in pep/*
do
    total_prot=`zcat $species | grep \> | wc -l`
    species=${species%.pep.all.fa.gz}
    species=${species#pep/}
    total_IDP=`expr $(cat output/${species}_IDP.csv | wc -l) - 1`
    species=`echo $species | awk '{match($0, /^_?[^_]+_[^_.]+/); print substr($0, RSTART, RLENGTH)}'` # Remove suffix
    echo $species, $total_prot, $total_IDP
done >> stat/IDP_percentage.csv
