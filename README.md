# slide-recon

the program use editdistance to find the UP

sample command to run count matrix:

python /path/to/slide-recon/run_bc_umi_extract.py -c cores -i input_directory -s sample_name

the code expects two fastq files with format sample_name_R[1-2]_001.fastq.gz inside input_direcory
the program starts with unzipping the fastqs inside input_dir then uses seqkit to split them in two multiple parts
number of parts is currently set equal to number of cores argument for processing, if there's not as many cores available to the program there will not be  but it can 
and process them separately 
