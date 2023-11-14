#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=80G
#SBATCH -t 1:55:00
#SBATCH -p short
#SBATCH -o recon_count_job_%A.out

echo 'indir =' $1
echo 'sample =' $2
echo 'max_anchors =' $3
echo 'max_targets =' $4

python ~/reconstruction/bc_umi_pipeline.py -c 16 -i $1 -s $2 -ma $3 -mt $4
#-l
