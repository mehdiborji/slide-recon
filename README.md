# slide-recon


## Sample command to run pipeline for generation of sparse count matrix from a pair of Illumina fastq file:

```
python /path/to/slide-recon/bc_umi_pipeline.py \
        -c cores \
        -i input_directory \
        -s sample_name \
        -ma max_anchors_bcs \
        -mt max_anchors_bcs \
        -r1 read1_structure \
        -r2 read2_structure
```

- This pipeline expects two fastq files with format `sample_name_R[1-2]_001.fastq.gz` inside `input_direcory`
The program starts with extracting and splitting the fastq files into a directory with path `input_direcory/sample_name/split`.
The number of reads per chunk can be modified with paramter within `bc_umi_utils.split_fastq_by_lines` function and has default of 10m reads.
The pipeline then uses `read_structure` information to process the splitted files in parallel and collect information counts of barcodes.
