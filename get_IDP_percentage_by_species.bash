#!/bin/bash

echo "FILENAME,TOTAL_IDP,NAME" > ./lib/IDP_percentage.csv

for file in ./pep/*
do
    acc=`basename $file`
    acc=${acc%.pep.all.fa.gz}
    total_IDP=`expr $(cat output/${acc}_IDP.csv | wc -l) - 1`
    species=`echo $acc | awk '{match($0, /^_?[^_]+_[^_.]+/); print substr($0, RSTART, RLENGTH)}'` # Remove suffix
    echo $acc,$total_IDP,$species
done >> ./lib/IDP_percentage.csv
