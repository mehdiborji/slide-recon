# slide-recon


## Sample command to run pipeline for generation of sparse count matrix from a pair of Illumina fastq files:

```
python /path/to/slide-recon/bc_umi_pipeline.py \
        -c cores \
        -i input_directory \
        -s sample_name \
        -ma max_r1_bcs \
        -mt max_r2_bcs \
        -r1 read1_structure \
        -r2 read2_structure
```

This pipeline expects two fastq files with format `sample_name_R[1-2]_001.fastq.gz` inside `input_direcory`. These are steps in the pipeline which uses functions within `bc_umi_utils.py` module.

1. Extracting and splitting the fastq files into an intermediate directory with path `input_direcory/sample_name/split`. The number of reads per chunk can be modified with paramter `lines` within `bc_umi_utils.split_fastq_by_lines` function and has default of 10m reads, or 40m lines.


2. Function `extract_bc_umi_dict` perform main pre-processing step on the raw reads
   - Fucntion `parse_read_struct` takes `read1_structure` and `read2_structure` from the input and converts the sequences into intervals for different landmark on the read, namely, barcode, UMI, polyT/A sequences, and other adapters placed in the bead oligo structure.
     
   - Function `edit_match` identifies reads containing appropriate landmarks by matching constant sequences of apadpters in known locations by calculating edit distances using `edlib` package.
     
   - Passing raw reads are written to the disk in `csv` format. Each line is a quadruple of R1_BC,R1_UMI,R2_BC,R2_UMI. also `json` files are stored for all R1_BC_R1_UMI combinations and R2_BC_R2_UMI combinations.
     
   - This is done in parallel for each fastq `part` using `multiprocessing` module.
  

3. Function `aggregate_dicts` aggregates R1 and R2 json files from all `parts` and stores in somewhat huge new json files.
   - This can be slow and final json file can be quite large due to potentially 100s of millions of raw barcodes. One solution is to ignore UMI information at this stage.
     
   - Function `aggregate_stat_dicts` aggregaitng barcode matching statistic for each of the adapters. This statitics is edit distance in the following format:
     
     `r1_adapter1_..._r1_adapterN__r2_adapter1_..._r2_adapterN`.

     For example -1_1__0_3 means r1_adapter1 didn't match, r1_adapter2 matched with 1 distance, and so on.
     
4. Function `whitelist_rankplot` uses R1 aggregated json files to estimate whitelist barcodes using a histogram-thresholding strategy on UMI counts of all raw barcode.
   - Same process is reapeated for R2 json files.
     
   - It also makes plots of the histrograms and rankplots underlying raw barcodes.
   
   - A dup-rate histogram for all whitelist barcodes is also generated. 

5. A second round of processing takes place on the quadruple csv files and they are converted into dictionary elements with key-value paring in 'R1_BC_whitelist':['R1_UMI','R2_BC_whitelist'] format.

6. Each of these `part` dictionaries is chunked into `batches` of `R1_BC_whitelist` barcodes, using function `save_barcode_batch_json`. Currently this is hardcoded to 30k barcodes in each `batch`.

7. Function `aggregate_barcode_batches` aggregates all parts for each of the batches and sotes them in new json files. 

8. Function `make_count_sparse_mtx_batch` takes each batch and stores it in a sparse count matrix.