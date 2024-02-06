#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=40G
#SBATCH -t 2:30:00
#SBATCH -p short
#SBATCH -o slide_gex_count_job_%A.out

echo 'indir =' $1
echo 'sample =' $2

python ~/reconstruction/gex_pipeline.py -c 20 -i $1 -s $2

