#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=48G
#SBATCH -t 0:50:00
#SBATCH -p short
#SBATCH -o slide_gex_count_job_%A.out

echo 'indir =' $1
echo 'sample =' $2

python ~/reconstruction/gex_pipeline.py -c 20 -i $1 -s $2

