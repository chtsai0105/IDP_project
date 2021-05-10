#!/bin/bash
#SBATCH -p batch
#SBATCH -J iupred2a_in_fungi
#SBATCH --time=3-00:00:00
#SABTCH -n 8
#SBATCH --mem=10g
#SBATCH --array=0-1013%8
#SBATCH --output=slurm_log/%x_%A_%a.out

files=( pep/* )
input=${files[$SLURM_ARRAY_TASK_ID]}
sample=${input%.pep.all.fa.gz}
sample=${sample#pep/}
output="output/IDR"

echo $SLURM_ARRAY_TASK_ID, $sample
./2_run_iupred_genome-wide.py --i $input --o $output
