import argparse
from multiprocessing import Pool
from matplotlib.backends.backend_pdf import PdfPages
import bc_umi_utils
import os
import edlib
import pandas as pd
import pysam
from tqdm import tqdm
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cores', type=str)
parser.add_argument('-i', '--indir', type=str)
parser.add_argument('-s', '--sample', type=str)
parser.add_argument('-l', '--limit', default=False, action='store_true')

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
limit = args.limit

######################################################

bc_umi_utils.split_fastq_by_lines(indir,sample,4e7)

######################################################
# the expected output is a 'one lane output from bcl-convert software'
######################################################

######################################################
# process the Demultiplex_Stats csv

demux = pd.read_csv(f'{indir}/Reports/Demultiplex_Stats.csv')

demux = demux[['SampleID', 'Index', '# Reads','% Reads', '# Perfect Index Reads',
       '# One Mismatch Index Reads', '# Two Mismatch Index Reads',
       '% Perfect Index Reads', '% One Mismatch Index Reads',
       '% Two Mismatch Index Reads']]
demux.columns=['id','index', 'nr', 'nrp', 'mm0','mm1','mm2','mm0p','mm1p','mm2p']
demux.nr = demux.nr/1e6
demux.mm0 = demux.mm0/1e6
demux.mm1 = demux.mm1/1e6
demux.mm2 = demux.mm2/1e6
demux = demux[~demux['index'].isna()]

######################################################
# process the Top_Unknown_Barcodes csv and select barcodes for further extraction

unkno = pd.read_csv(f'{indir}/Reports/Top_Unknown_Barcodes.csv')
unkno = unkno[['index', 'index2', '# Reads', '% of Unknown Barcodes']].copy()
unkno.columns=['index', 'index2', 'nr', '%']
unkno.nr = unkno.nr/1e6
unkno['merge_idx'] = unkno['index'] + '+' + unkno['index2']

# we have two 8nt indices we have observed to contain many lost reads

target_bcs = ['CGTACTAG', 'CTCTCTAC'] 

bc_sample_dict = {}
for bc in target_bcs:
    
    bc_sample_dict[bc] = demux[demux['index'].str.contains(bc)]['id'].tolist()[0]

######################################################
# calculat distance of Top_Unknown_Barcodes

bc = 'AGATCTCGGT' # generic sequence expected at place of index2 for single index libraries
eds = []
                    
for err_bc in unkno['index2']:
    ed = edlib.align(bc,err_bc)['editDistance']
    eds.append(ed)
unkno[f'{bc}_dist'] = eds

for bc in target_bcs:
    eds = []
    for err_bc in unkno['index']:
        ed = edlib.align(bc,err_bc[:-2])['editDistance']
        eds.append(ed)
    unkno[f'{bc}_dist']  = eds

######################################################
# select Unknown_Barcodes with small distance to index2 and each of target_bcs for further extraction
                    
accept_bcs = {}

sub_unkno = unkno[ (unkno.AGATCTCGGT_dist<=3) & (unkno.CTCTCTAC_dist<=4) & (unkno.CGTACTAG_dist>4)]
print(sub_unkno['nr'].sum())
accept_bcs['CTCTCTAC'] = sub_unkno['merge_idx'].tolist()

sub_unkno = unkno[ (unkno.AGATCTCGGT_dist<=3) & (unkno.CTCTCTAC_dist>4) & (unkno.CGTACTAG_dist<=4)]
print(sub_unkno['nr'].sum())
accept_bcs['CGTACTAG'] = sub_unkno['merge_idx'].tolist()

######################################################

parts = bc_umi_utils.find_sub_fastq_parts(indir,sample)
args = [(indir,sample,part,limit) for part in parts]

N_read_extract = 10000

def write_fastq_pair_rescue(R1_rescue, R2_rescue, r1, r2):
    
    R1_rescue.write(f'@{r1.name}\n')
    R1_rescue.write(f'{r1.sequence}\n')
    R1_rescue.write('+\n')
    R1_rescue.write(f'{r1.quality}\n')

    R2_rescue.write(f'@{r2.name}\n')
    R2_rescue.write(f'{r2.sequence}\n')
    R2_rescue.write('+\n')
    R2_rescue.write(f'{r2.quality}\n')

def extract_rescue(indir,sample,part,limit):
    
    i = 0
    
    R1_fastq = f'{indir}/{sample}/split/{sample}_R1.part_{part}.fastq'
    R2_fastq = f'{indir}/{sample}/split/{sample}_R2.part_{part}.fastq'
    
    target_fastqs = {}
    
    for bc in target_bcs:
        target_fastqs[f'R1_{bc}'] = R1_fastq.replace(f'split/{sample}',bc_sample_dict[bc] + '_rescue')
        target_fastqs[f'R2_{bc}'] = R2_fastq.replace(f'split/{sample}',bc_sample_dict[bc] + '_rescue')
    
    if os.path.isfile(target_fastqs[f'R1_{bc}']):
        print(target_fastqs[f'R1_{bc}'],' exists, skip')
    else:
        print(target_fastqs[f'R1_{bc}'],' does not exist, will extract')
        
        target_ios = {}
        
        for bc in target_bcs:
            target_ios[f'R1_{bc}'] = open(R1_fastq.replace(f'split/{sample}',bc_sample_dict[bc] + '_rescue'), 'w')
            target_ios[f'R2_{bc}'] = open(R2_fastq.replace(f'split/{sample}',bc_sample_dict[bc] + '_rescue'), 'w')

        with pysam.FastxFile(R1_fastq) as R1, pysam.FastxFile(R2_fastq) as R2:
            for r1, r2 in tqdm(zip(R1, R2)):
                i += 1;
                idx_pair = r1.comment.split('N:0:')[1]
                for bc in target_bcs:
                    if idx_pair in accept_bcs[bc]:

                        write_fastq_pair_rescue(target_ios[f'R1_{bc}'],target_ios[f'R2_{bc}'],r1,r2)

                if i>N_read_extract and limit: break

        for bc in target_bcs:
            target_ios[f'R1_{bc}'].close()
            target_ios[f'R2_{bc}'].close()
    
    for file in list(target_fastqs.values()):
        
        file_gz = file.replace('.fastq','.fastq.gz')
        
        if os.path.isfile(file_gz):
            print(file_gz,' exists, skip')
        else:
            print(file_gz,' does not exist, will zip')
            subprocess.call(f'pigz -c {file} > {file_gz}',shell=True)
            

pool = Pool(int(cores))
results = pool.starmap(extract_rescue, args)
pool.close()
pool.join()

for bc in target_bcs:
    subprocess.call(f'cat {indir}/{sample}/{bc_sample_dict[bc]}_rescue_R1.part*.fastq.gz > {indir}/{bc_sample_dict[bc]}_rescue_R1_001.fastq.gz',shell=True)
    subprocess.call(f'cat {indir}/{sample}/{bc_sample_dict[bc]}_rescue_R2.part*.fastq.gz > {indir}/{bc_sample_dict[bc]}_rescue_R2_001.fastq.gz',shell=True)