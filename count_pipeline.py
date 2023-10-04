import pysam
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import edlib
import argparse
import csv
from pathlib import Path
import json
import subprocess
from multiprocessing import Pool
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
from anndata import AnnData
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument('--cores', type=str)
parser.add_argument('--indir', type=str)
parser.add_argument('--sample', type=str)
parser.add_argument('--limit', default=False, action='store_true')
#parser.add_argument('--outdir', type=str)
#parser.add_argument('--barcodes', type=str)
#parser.add_argument('--split', default=False, action='store_true')
#parser.add_argument('--mode', type=str)

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
limit = args.limit

#outdir = args.outdir
#barcodes = args.barcodes
#split = args.split
#mode = args.mode

UP_seq='TCTTCAGCGTTCCCGAGA'
ad_seq=[('N', 'A'), ('N', 'T'), ('N', 'G'), ('N', 'C')]
#print(UP_seq)

"""
fastq_dir='/Users/mborji/reconstruct/230803/fastq/'

fastq_dir='/Users/mborji/reconstruct/all_fastq/'

files=sorted(os.listdir(fastq_dir))

samples=np.unique([f.split('_R')[0] for f in files if 'fastq.gz' in f])

sample=samples[0]

"""

N_read_extract=100000

print(N_read_extract)

def split_fastq(indir,sample,cores):
    
    splitted_file = f'{indir}/{sample}/split/{sample}_R1_001.part_001.fastq.gz'
    
    if os.path.isfile(splitted_file):
        print(splitted_file,' splitted fastq exists, skip splitting')
    else:
        print(splitted_file,' splitted fastq does not exist')
        R1 = f'{indir}/{sample}_R1_001.fastq.gz'
        R2 = f'{indir}/{sample}_R2_001.fastq.gz'
        subprocess.call(['seqkit', 'split2' , '-1', R1, '-2', R2, '-p', str(cores), '-f', '-O', f'{indir}/{sample}/split/'])

def unzip_split_fastq(indir,sample,cores):
    
    splitted_file = f'{indir}/{sample}/split/{sample}_R1_001.part_001.fastq'
    
    if os.path.isfile(splitted_file):
        print(splitted_file,' splitted fastq exists, skip splitting')
    else:
        print(splitted_file,' splitted fastq does not exist')
        R1 = f'{indir}/{sample}_R1_001.fastq.gz'
        R2 = f'{indir}/{sample}_R2_001.fastq.gz'
        
        unR1 = R1.replace('.gz','')
        unR2 = R2.replace('.gz','')
        
        command=f'zcat {R1} > {unR1}'
        subprocess.call(command,shell=True)
        command=f'zcat {R2} > {unR2}'
        subprocess.call(command,shell=True)

        subprocess.call(['seqkit', 'split2' , '-1', unR1, '-2', unR2, '-p', 
                         str(cores), '-f', '-O', f'{indir}/{sample}/split/'])
        
        command=f'rm {unR1} {unR2}'
        subprocess.call(command,shell=True)
        
        #subprocess.call(['pigz', '-f', f'{indir}/{sample}/split/*.fastq'])

        
def seq_counter(seq_dict,seq_instance):
    
    if seq_dict.get(seq_instance) is None:
        seq_dict[seq_instance] = 1
    else:
        seq_dict[seq_instance]+= 1

def quad_dict_store(quad_dict,quad_key,quad_items):
    
    if quad_dict.get(quad_key) is None:
        quad_dict[quad_key] = [quad_items]
    else:
        quad_dict[quad_key].extend([quad_items])
        
def UP_edit_pass(read_seq,max_dist):
    edit=edlib.align(read_seq[8:26],UP_seq,'HW','path',max_dist,ad_seq)
    boolean_pass = edit['editDistance']>=0 and edit['editDistance']<=max_dist
    return(boolean_pass)

def seq_slice(read_seq):
    bc = read_seq[:8]+read_seq[26:32]
    umi = read_seq[33:41]
    return(bc,umi)

def outfile_from_in(infile,suffix):
    split_root=os.path.dirname(infile)
    split_part=os.path.basename(infile).split('.')[1]
    outfile=f'{split_root}/{sample}.{split_part}_{suffix}.json'
    return(split_root,split_part,outfile)

def extract_bc_umi_dict(R1_fastq,R2_fastq):
    
    i=0;max_dist=2

    split_root, split_part, targets_json=outfile_from_in(R1_fastq,'targets')
    split_root, split_part, anchors_json=outfile_from_in(R1_fastq,'anchors')
    
    if os.path.isfile(anchors_json):
        print(anchors_json,' exists, skip')
        return
    #if os.path.isfile(anchors_json):
    #    print(anchors_json,' exists, skip extracting')
    #    return
    
    targets_dict = {}
    anchors_dict = {}
    
    with pysam.FastxFile(R1_fastq) as R1, pysam.FastxFile(R2_fastq) as R2:
        for r1, r2 in tqdm(zip(R1, R2)):
            i+=1
            seq1 = r1.sequence
            seq2 = r2.sequence
            
            if UP_edit_pass(seq1,max_dist) and UP_edit_pass(seq2,max_dist):
                
                a_bc,a_umi=seq_slice(seq1)
                t_bc,t_umi=seq_slice(seq2)
                quad_dict_store(targets_dict,t_bc,t_umi)
                quad_dict_store(anchors_dict,a_bc,a_umi)
                
            if i>N_read_extract and limit: break
            
    with open(targets_json, 'w') as json_file:
        json.dump(targets_dict, json_file)
    with open(anchors_json, 'w') as json_file:
        json.dump(anchors_dict, json_file)

def aggregate_dicts(indir,sample,position): 
    
    dir_split=f'{indir}/{sample}/split/'
    files=os.listdir(dir_split)
    jsons = sorted([f for f in files if f'{position}.json' in f])
    
    agg_read_csv=f'{indir}/{sample}/{sample}_agg_read_cnt_{position}.csv'
    
    if os.path.isfile(agg_read_csv):
        print(agg_read_csv,' exists, skip')
        return
    
    data_agg={}
    for i in tqdm(range(len(jsons))):
        with open(f'{dir_split}{jsons[i]}', 'r') as json_file:
            data_sub = json.load(json_file)
            print(jsons[i],len(data_sub))
            for k in data_sub:
                if data_agg.get(k) is not None:
                    data_agg[k].extend(data_sub[k])
                else:
                    data_agg[k]=data_sub[k]
    read_dict={}
    umi_dict={}
    total_reads=0
    for k in tqdm(data_agg):
        reads=len(data_agg[k])
        total_reads+=reads
        if reads>=5:
            read_dict[k]=reads
            umi_dict[k]=len(set(data_agg[k]))
            
    umi_cnt=pd.Series(umi_dict)
    read_cnt=pd.Series(read_dict)
    
    read_cnt.to_csv(agg_read_csv)
    umi_cnt.to_csv(agg_read_csv.replace('read','umi'))
    
def whitelist_rankplot(indir,sample,position):
    
    #dup_rate_file=f'{indir}/{sample}/{sample}_{position}_duprate.pdf'
    #if os.path.isfile(dup_rate_file):
    #    print(dup_rate_file,' exists, skip')
    #    return
    
    read_cnt=pd.read_csv(f'{indir}/{sample}/{sample}_agg_read_cnt_{position}.csv')
    umi_cnt=pd.read_csv(f'{indir}/{sample}/{sample}_agg_umi_cnt_{position}.csv')
    umi_cnt.columns=['bc','umi_cnt']
    read_cnt.columns=['bc','read_cnt']
    agg_bcs=pd.merge(umi_cnt,read_cnt,left_on='bc',right_on='bc',how='inner')
    agg_bcs['log10_read_cnt']=np.log10(agg_bcs['read_cnt'])
    agg_bcs['log10_umi_cnt']=np.log10(agg_bcs['umi_cnt'])
    agg_bcs['dup_rate']=agg_bcs['read_cnt']/agg_bcs['umi_cnt']
    agg_bcs=agg_bcs.sort_values(by='umi_cnt',ascending=False)
    
    sub=agg_bcs.iloc[20:100000].copy()  # select top 100k bc except first 20
    x = np.histogram(sub.log10_umi_cnt, 100) # fit a histogram
    smooth = gaussian_filter1d(x[0], 3) # smooth histogram
    peak_idx,_=find_peaks(-smooth) # find the local minimum
    print(peak_idx,x[1][:-1][peak_idx])
    mean_hist=(x[1][1:][peak_idx]+x[1][:-1][peak_idx])/2 # take the mid point of point before and after
    
    mean_hist=mean_hist[-1] # take the last value in list of local minima (could be more than one)

    wl_df=agg_bcs[agg_bcs.log10_umi_cnt>=mean_hist].copy()
    wl_df.to_csv(f'{indir}/{sample}/{sample}_{position}_wl.csv.gz',compression='infer')
    #wl_reads=wl_df.read_cnt.sum()
    white_list_size=wl_df.shape[0]
    
    plt.figure(figsize=(4,3))
    log10_ranks=np.log10(np.arange(1,len(agg_bcs)+1))
    log10_cnts=agg_bcs.log10_umi_cnt
    plt.plot(log10_ranks,log10_cnts)#,label='Rank Plot of Reads')
    plt.xlabel('Log10 Ranks')
    plt.ylabel('Log10 UMI Counts')
    plt.title(f'{sample} {position}\n {white_list_size} white listed')
    plt.plot([0, log10_ranks[-1]], [mean_hist, mean_hist], linewidth=1,label='log10 threshold',c='tab:green')
    log10_wl=np.log10(white_list_size)
    plt.plot([log10_wl, log10_wl], [log10_cnts.min(), log10_cnts.max()], linewidth=1,label='log10 size',c='tab:orange')
    plt.legend(loc="best");
    
    qc_pdfs.savefig(bbox_inches='tight')
    #plt.savefig(f'{indir}/{sample}/{sample}_{position}_rankplot.pdf',bbox_inches='tight');
    
    plt.figure(figsize=(4,3))
    plt.plot(x[1][:-1],x[0], label='Raw Histogram')
    plt.plot(x[1][:-1],smooth, label='Gaussian Smoothed')
    plt.xlabel('Log10 UMI Counts')
    plt.ylabel('Bin Height')
    plt.title(f'{sample} {position}')
    plt.plot([mean_hist, mean_hist], [0, np.max(x[0])], linewidth=2,label='Whitelist Threshold')
    plt.legend(loc="best");
    
    qc_pdfs.savefig(bbox_inches='tight')
    #plt.savefig(f'{indir}/{sample}/{sample}_{position}_histogram.pdf',bbox_inches='tight');
    
    plt.figure(figsize=(2,2))
    mean=wl_df.dup_rate.mean()
    n_std=3
    width=wl_df.dup_rate.std()*n_std
    ticks=np.linspace(mean-width,mean+width,n_std)
    sns.histplot(wl_df[(wl_df.dup_rate>mean-width) & (wl_df.dup_rate<mean+width) ].dup_rate,bins=50)
    plt.xticks(ticks)
    
    qc_pdfs.savefig(bbox_inches='tight')
    #plt.savefig(f'{indir}/{sample}/{sample}_{position}_duprate.pdf',bbox_inches='tight');

def extract_quad_dict(indir,sample,part):
    
    i = 0;max_dist = 2;quad_dict = {}
    
    part = f'part_{str(part).zfill(3)}'
    
    R1_fastq = f'{indir}/{sample}/split/{sample}_R1_001.{part}.fastq'
    R2_fastq = f'{indir}/{sample}/split/{sample}_R2_001.{part}.fastq'
    
    quads_json = f'{indir}/{sample}/split/{sample}.{part}_quads.json'
    if os.path.isfile(quads_json):
        print(quads_json,' exists, skip')
        return
    
    a_white = pd.read_csv(f'{indir}/{sample}/{sample}_anchors_wl.csv.gz')['bc']
    t_white = pd.read_csv(f'{indir}/{sample}/{sample}_targets_wl.csv.gz')['bc']

    a_dict = {}
    for bc in a_white: a_dict[bc] = []
    t_dict = {}
    for bc in t_white: t_dict[bc] = []
    
    with pysam.FastxFile(R1_fastq) as R1, pysam.FastxFile(R2_fastq) as R2:
        for r1, r2 in tqdm(zip(R1, R2)):
            i+=1
            
            seq1 = r1.sequence
            seq2 = r2.sequence
            #if UP_edit_pass(seq1,max_dist) and UP_edit_pass(seq2,max_dist):
            a_bc,a_umi=seq_slice(seq1)
            t_bc,t_umi=seq_slice(seq2)

            if (a_dict.get(a_bc) is not None) and (t_dict.get(t_bc) is not None):
                quad_dict_store(quad_dict,a_bc,[a_umi,t_bc])
            if i>N_read_extract and limit: break
            
    with open(quads_json, 'w') as json_file:
        json.dump(quad_dict, json_file)
        
def find_sub_fastq_pairs(indir,sample):

    R1s=sorted([f for f in os.listdir(f'{indir}/{sample}/split/') if '_R1_001.part' in f])
    R2s=[f.replace('_R1_','_R2_') for f in R1s]
    pairs=[]

    for i in range(len(R1s)): 
        pairs.append((f'{indir}/{sample}/split/{R1s[i]}',f'{indir}/{sample}/split/{R2s[i]}'))
        
    return pairs

def make_count_mtx(indir,sample):
    
    position='quads'

    dir_split=f'{indir}/{sample}/split/'
    files=os.listdir(dir_split)
    jsons = sorted([f for f in files if f'{position}.json' in f])

    data_agg={}

    for i in tqdm(range(len(jsons))):
        with open(f'{dir_split}{jsons[i]}', 'r') as json_file:
            data_sub = json.load(json_file)
            print(jsons[i],len(data_sub))
            for k in data_sub:
                if data_agg.get(k) is not None:
                    data_agg[k].extend(data_sub[k])
                else:
                    data_agg[k]=data_sub[k]

    a_white = pd.read_csv(f'{indir}/{sample}/{sample}_anchors_wl.csv.gz',index_col=1)#['bc']
    t_white = pd.read_csv(f'{indir}/{sample}/{sample}_targets_wl.csv.gz',index_col=1)#['bc']
    
    counts_np = np.zeros( (len(a_white),len(t_white)) )
    counts_df = pd.DataFrame(counts_np, index=a_white.index, columns=t_white.index)

    all_list=[]
    for a_bc in tqdm(a_white.index):

        umi_tbc=data_agg[a_bc]
        umi_bc_dic={}
        for a in umi_tbc:
            if umi_bc_dic.get(a[0]) is not None:
                umi_bc_dic[a[0]].append(a[1])
            else:
                umi_bc_dic[a[0]]=[a[1]]

        t_bc_cnt={}
        for k in umi_bc_dic:
            umi_reads=len(umi_bc_dic[k])
            if umi_reads>2:
                uni_t_bc=set(umi_bc_dic[k])
                if len(uni_t_bc)>1:
                    bcs,cnts=np.unique(umi_bc_dic[k],return_counts=True)
                    if np.max(cnts/umi_reads)>.74:
                        t_bc=bcs[np.argmax(cnts/umi_reads)]
                        seq_counter(t_bc_cnt,t_bc)
                else:
                    t_bc=list(uni_t_bc)[0]
                    seq_counter(t_bc_cnt,t_bc)

        counts_df.loc[a_bc,list(t_bc_cnt.keys())]=list(t_bc_cnt.values())
        
    counts_df=AnnData(counts_df,dtype='float32')
    counts_df.X = csr_matrix(counts_df.X)
    counts_df.write_h5ad(f'{indir}/{sample}/{sample}_counts_filtered.h5ad')

#if __name__ == '__main__':

unzip_split_fastq(indir,sample,cores)

#split_fastq(indir,sample,cores)

######################################################
args= find_sub_fastq_pairs(indir,sample)
[print(a) for a in args]
pool = Pool(int(cores))
results = pool.starmap(extract_bc_umi_dict, args)
pool.close()
pool.join()
######################################################
aggregate_dicts(indir,sample,'anchors')
aggregate_dicts(indir,sample,'targets')

qc_pdfs = PdfPages(f'{indir}/{sample}/{sample}_QC.pdf')
whitelist_rankplot(indir,sample,'anchors')
whitelist_rankplot(indir,sample,'targets')
qc_pdfs.close()
######################################################
args=[(indir,sample,i) for i in range(1,int(cores)+1)]
[print(a) for a in args]
pool = Pool(int(cores))
results = pool.starmap(extract_quad_dict, args)
pool.close()
pool.join()
######################################################

make_count_mtx(indir,sample)
#a_white = pd.read_csv(f'{indir}/{sample}/{sample}_anchors_wl.csv.gz',index_col=1)#['bc']

#t_white = pd.read_csv(f'{indir}/{sample}/{sample}_targets_wl.csv.gz',index_col=1)#['bc']

"""
targets_json=outfile_from_in(args[0][0],'targets')

if os.path.isfile(targets_json):
    print(targets_json,' targets json exist, skip extracting')
else:
    print(targets_json,' splitted fastq does not exist')
    

pool = Pool(int(cores/2))
#results = pool.starmap(extract_bc_umi, args)
results = pool.starmap(extract_bc_umi_dict, args)
pool.close()
pool.join()

"""