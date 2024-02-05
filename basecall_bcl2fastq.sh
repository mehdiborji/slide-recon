#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 0:10:00
#SBATCH -p short
#SBATCH -o basecall_job_%A.out

module load bcl2fastq/2.20.0.422

echo 'input_folder =' $1
echo 'output_folder =' $2
echo 'sheet_file =' $3

bcl2fastq -R $1 -o $2 --sample-sheet $3