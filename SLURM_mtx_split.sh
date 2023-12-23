#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 2:00:00
#SBATCH -p short
#SBATCH -o mtx_job_%A.out

echo 'indir =' $1

python ~/reconstruction/mtx_dropseq_split.py --c 20 --i $1