#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 18:00:00
#SBATCH -p priority
#SBATCH -o umap_job_%A.out

echo 'indir =' $1
echo 'sample =' $2
echo 'adata_name =' $3
echo 'metric =' $4
echo 'n_neighbors =' $5
echo 'min_dist =' $6
#echo 'subset =' $3
#echo 'threshold =' $4

python ~/reconstruction/umap_reduction.py --c 16 --i $1 --sample $2 --adata_name $3 --metric $4 --n_neighbors $5 --min_dist $6
#python ~/reconstruction/umap_reduction.py --c 16 --i $1 --sample $2 --adata_name $3
#python ~/reconstruction/umap_reduction.py --c 16 --i $1 --sample $2 --subset $3 --threshold $4
