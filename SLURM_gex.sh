#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 0:20:00
#SBATCH -p short
#SBATCH -o slide_recon_count_job_%A.out

echo 'indir =' $1
echo 'sample =' $2

python ~/reconstruction/gex_pipeline.py -c 20 -i $1 -s $2

