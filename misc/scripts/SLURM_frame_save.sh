#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 0:10:00
#SBATCH -p priority
#SBATCH -o umap_job_%A.out

python ~/reconstruction/frame_save.py --c 20
