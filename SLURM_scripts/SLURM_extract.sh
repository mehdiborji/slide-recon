#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH -t 0:50:00
#SBATCH -p priority
#SBATCH -o slide_recon_count_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'max_anchors =' $3
echo 'max_targets =' $4

python ~/reconstruction/bc_umi_pipeline.py -c 4 -i $1 -s $2 -ma $3 -mt $4 -r1 $5 -r2 $6
#--limit
