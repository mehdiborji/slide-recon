#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=2G
#SBATCH -t 0:50:00
#SBATCH -p short
#SBATCH -o nova_rescue_job_%A.out

echo 'indir =' $1
echo 'sample =' $2

#./SLURM_rescue.sh /n/scratch/users/m/meb521/Nova_L3_final Undetermined_S0_L003
python ~/reconstruction/rescue_pipeline.py -c 20 -i $1 -s $2

