#!/bin/bash

#$ -N slide_recon_count.out
#$ -l os=RedHat7
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -pe smp 4
#$ -binding linear:4
#$ -terse
#$ -notify
#$ -R y
#$ -j y
#$ -m eas

echo 'indir =' $1
echo 'sample =' $2
echo 'max_anchors =' $3
echo 'max_targets =' $4

python ~/slide-recon/bc_umi_pipeline.py -c 8 -i $1 -s $2 -ma $3 -mt $4 -r1 $5 -r2 $6