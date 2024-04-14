#!/bin/bash

#$ -N slide_recon_count.out
#$ -l os=RedHat7
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -binding linear:1
#$ -terse
#$ -notify
#$ -R y
#$ -j y
#$ -m eas

echo 'indir =' $1
echo 'sample =' $2
echo 'adata_name =' $3
echo 'metric =' $4
echo 'n_neighbors =' $5
echo 'min_dist =' $6
echo 'spread =' $7

source /broad/software/scripts/useuse
reuse Anaconda3
#reuse -q Anaconda3
#source activate /home/unix/chu/anaconda3/envs/slidelock


python ~/slide-recon/umap_reduction.py --c 8 --i $1 --sample $2 --adata_name $3 \
--metric $4 --n_neighbors $5 --min_dist $6  --spread $7