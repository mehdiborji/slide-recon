# slide-recon

the program use editdistance to find the UP

sample command to run count matrix:

python /path/to/slide-recon/run_bc_umi_extract.py -c cores -i input_directory -s sample_name

sbatch ~/reconstruction/SLURM_extract.sh . H3_1_merge 400000 90000

the code expects two fastq files with format sample_name_R[1-2]_001.fastq.gz inside input_direcory
the program starts with unzipping the fastqs inside input_dir then uses seqkit to split them in two multiple parts
number of parts is currently set equal to number of cores argument for processing, if there's not as many cores available to the program there will not be  but it can 
and process them separately 


## fiducial_seq_read.py
match barcode from fastq to in situ sequencing whitelist, extract barcode and output anchor-target-umi pair list

run by
~ python fiducial_seq_read.py  -d <231104> -s <c58_16> -p <Puck_230818_16>


## plot_groundtruth.py
plot diffusion property based on ground truth
input matched reads from fiducial_seq_read.py

1. plot distribution of umi count and unique diffused bc per anchor and target bead (func: plot_cnt_distribution)
2. plot diffusion pattern of anchor and target on a 2D space, save as pdf (func: plot_diffusion)
3. plot diffusion on x and y axis separately (func: plot_1d_diffusion)
4. average over diffusion and kit kde on x and y axis (func: plot_kde_average)

output png and pdf in sample folder

run by
~ python plot_groundtruth.py  -d <231104> -s <c58_16> -a <V9A30> -t <V10>


## fiducial_seq_blind_whitelist.py
extracting barcode from fastq
output anchor-target-umi pair list

1. extract barcode by parsing fastq, check UP first for distance <=3 (func: barcode_extract)
2. plot barcode rank for both anchor and target
3. set threshold as min reads per barcode for a smaller barcode collections, use UMICluster ('directional') from UMI-cools to collapse and the collapsed result serves as barcode whitelist. Also match barcodes below the threshold to the whitelist as collapsing. Generate dict with anchor-target-umi pair (func: bc_collapsing, umi_collapsing)
4. wirting anchor-target-umi pair to csv.gz and also report parsing stats (func: write_blind)

run by
~ python fiducial_seq_blind_whitelist.py  -d <231104> -s <H2_15>