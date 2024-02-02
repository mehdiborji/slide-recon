#!/bin/bash

#SBATCH -c 20
#SBATCH --mem=80G
#SBATCH -t 0:25:00
#SBATCH -p short
#SBATCH -o alignment_job_%A.out

genome_dir=$1
out_name=$2
input_fastq_r1=$3
input_fastq_r2=$4
input_whitelist=$5

echo $input_fastq_r1
echo $input_fastq_r2
echo $genome_dir
echo $out_name
echo $input_whitelist

STAR \
--runThreadN 20 \
--readFilesIn $input_fastq_r1 $input_fastq_r2 \
--genomeDir $genome_dir \
--outFileNamePrefix $out_name \
--soloCBwhitelist $input_whitelist \
--soloType CB_UMI_Simple \
--soloCBstart 1 \
--soloCBlen 15 \
--soloUMIstart 16 \
--soloUMIlen 9 \
--soloFeatures Gene Velocyto GeneFull \
--soloUMIdedup Exact \
--soloMultiMappers EM \
--soloCellFilter None \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes AS CR UR CB UB GX GN \
--soloOutFormatFeaturesGeneField3 - \
--soloCellReadStats Standard
#--readFilesCommand zcat
#--alignIntronMax 1 \
#--readMapNumber $5 \
#--outSAMmode NoQS \
#--outSAMattributes AS \
#--outFilterMultimapNmax 10 \
#samtools sort -@16 -o "$out_name".bam "$out_name"Aligned.out.sam
#samtools index -@16 "$out_name".bam
#samtools view -h -b -F 2308 -@16 "$out_name".bam > "$out_name"_pri.bam
#samtools index -@16 "$out_name"_pri.bam

#rm "$out_name"*out*

