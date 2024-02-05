#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=40G
#SBATCH -t 0:40:00
#SBATCH -p priority
#SBATCH -o basecall_job_%A.out

#source /broad/software/scripts/useuse
#use .bcl2fastq2-v2.20.0

module load bclConvert/4.2.7

echo 'input_folder =' $1
echo 'output_folder =' $2
echo 'sheet_file =' $3

bcl-convert \
--bcl-input-directory $1 \
--output-directory $2 \
--sample-sheet $3 \
--force
#--bcl-num-conversion-threads 20
#--bcl-num-compression-threads 2 \
#--bcl-num-decompression-threads 4
