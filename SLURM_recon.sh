#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=6G
#SBATCH -t 0:10:00
#SBATCH -p priority
#SBATCH -o recon_job_%A.out

echo 'indir =' $1
echo 'sample =' $2

python ~/reconstruction/run_bc_umi_extract.py --c 16 --i $1 --s $2
#--l
