#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH -t 0:40:00
#SBATCH -p short
#SBATCH -o umap_job_%A.out

echo 'indir =' $1
echo 'sample =' $2
echo 'adata_name =' $2
#echo 'subset =' $3
#echo 'threshold =' $4

python ~/reconstruction/umap_reduction.py --c 16 --i $1 --sample $2 --adata_name $3
#python ~/reconstruction/umap_reduction.py --c 16 --i $1 --sample $2 --subset $3 --threshold $4
