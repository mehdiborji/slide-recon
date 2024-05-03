#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=80G
#SBATCH -t 30:00:00
#SBATCH -p medium
#SBATCH -o umap_job_%A.out
#SBATCH --account=chen_fec176

echo 'indir =' $1
echo 'sample =' $2
echo 'adata_name =' $3
echo 'metric =' $4
echo 'n_neighbors =' $5
echo 'min_dist =' $6
echo 'spread =' $7
#echo 'subset =' $3
#echo 'threshold =' $4

python ~/reconstruction/umap_reduction.py --c 20 --i $1 --sample $2 --adata_name $3 \
--metric $4 --n_neighbors $5 --min_dist $6  --spread $7
