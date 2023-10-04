# slide-recon

sample command to run is
python /path/to/slide-recon/run_bc_umi_extract.py -c cores -i input_directory -s sample_name
the code expects two fastq files with format sample_name_R[1-2]_001.fastq.gz inside input_direcory
the program starts with unzipping the fastqs inside input_dir
then uses sektk to split them
