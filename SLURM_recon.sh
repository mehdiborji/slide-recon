#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 1:10:00
#SBATCH -p short
#SBATCH -o recon_count_job_%A.out

echo 'indir =' $1
echo 'sample =' $2

python ~/reconstruction/run_bc_umi_extract.py --c 16 --i $1 --s $2
#--l
