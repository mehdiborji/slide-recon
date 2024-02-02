import pysam
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import edlib
import csv
from pathlib import Path
import json
import subprocess
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
from anndata import AnnData
import re
from collections import defaultdict

UP_seq = 'TCTTCAGCGTTCCCGAGA'

filter_poly = True

N_read_extract = 500000

print(N_read_extract)

def split_fastq_by_lines(indir,sample,lines=4e6):
    
    splitted_file = f'{indir}/{sample}/split/{sample}_R1.part_000.fastq'
    
    if os.path.isfile(splitted_file):
        print(splitted_file,' splitted fastq exists, skip splitting')
    else:
        print(splitted_file,' splitted fastq does not exist')
        
        R1 = f'{indir}/{sample}_R1_001.fastq.gz'
        R2 = f'{indir}/{sample}_R2_001.fastq.gz'
        
        split_dir = f'{indir}/{sample}/split'
        if not os.path.exists(split_dir):
            os.makedirs(split_dir)
            print(f'{split_dir} created')
        else:
            print(f'{split_dir} already exists')
        
        split_R1_name = f'{split_dir}/{sample}_R1.part_'
        split_R2_name = f'{split_dir}/{sample}_R2.part_'
        
        # zcat and split by n=lines (4x number of reads) add a suffix with 3 digits and prefix of 'split_R1_name'
        command_R1 = f'zcat {R1} | split -a 3 -l {int(lines)} -d --additional-suffix=.fastq - {split_R1_name}'
        command_R2 = command_R1.replace('_R1','_R2')

        subprocess.call(f'{command_R1} & {command_R2}', shell=True)

def seq_counter(seq_dict,seq_instance):
    
    if seq_dict.get(seq_instance) is None:
        seq_dict[seq_instance] = 1
    else:
        seq_dict[seq_instance] += 1

def quad_dict_store(quad_dict,quad_key,quad_items):
    
    if quad_dict.get(quad_key) is None:
        quad_dict[quad_key] = [quad_items]
    else:
        quad_dict[quad_key].extend([quad_items])
        
def UP_edit_pass(read_seq,max_dist):
    
    edit = edlib.align(read_seq[8:26],UP_seq,'HW','distance',max_dist)
    ed_dist = edit['editDistance']
    boolean_pass = ed_dist >= 0 and ed_dist <= max_dist
    
    return(boolean_pass, ed_dist)

def seq_slice(read_seq):
    
    bc = read_seq[:8]+read_seq[26:33]
    umi = read_seq[33:42]
    return(bc, umi)

def find_sub_fastq_parts(indir,sample):
    
    pattern = re.compile(r'_R1.part_(.*?)\.fastq')
    all_files = os.listdir(f'{indir}/{sample}/split/')
    parts = sorted([f.split('.part_')[1].split('.fastq')[0] for f in all_files if pattern.search(f)])
    parts = sorted(np.unique([f.split('.part_')[1][:3] for f in all_files if pattern.search(f)])) # part + 3 digits because we did split suffix with 3 digits
    
    return parts

def extract_bc_umi_dict(indir,sample,part,limit):
    
    i = 0
    max_dist = 2
    
    R1_fastq = f'{indir}/{sample}/split/{sample}_R1.part_{part}.fastq'
    R2_fastq = f'{indir}/{sample}/split/{sample}_R2.part_{part}.fastq'
    
    anchors_json = f'{indir}/{sample}/split/{sample}.part_{part}_anchors.json'
    targets_json = f'{indir}/{sample}/split/{sample}.part_{part}_targets.json'
    
    anchors_umi_len_json = f'{indir}/{sample}/split/{sample}.part_{part}_anchors_umi_len.json'
    targets_umi_len_json = f'{indir}/{sample}/split/{sample}.part_{part}_targets_umi_len.json'
    
    anchor_edits_json = f'{indir}/{sample}/split/{sample}.part_{part}_anchor_edits.json'
    target_edits_json = f'{indir}/{sample}/split/{sample}.part_{part}_target_edits.json'
    
    if os.path.isfile(anchors_json):
        print(anchors_json,' exists, skip')
        return
    
    anchors_dict = {}
    targets_dict = {}
    
    anchor_edits_dict = {}
    target_edits_dict = {}
    
    anchors_umi_len_dict = {}
    targets_umi_len_dict = {}
    
    with pysam.FastxFile(R1_fastq) as R1, pysam.FastxFile(R2_fastq) as R2:
        for r1, r2 in tqdm(zip(R1, R2)):
            
            i+=1
            
            seq1 = r1.sequence
            seq2 = r2.sequence
            
            len1 = len(seq1)
            len2 = len(seq2)
            
            if filter_poly:
                
                polyT_cnt = seq1[42:50].count('T')
                polyA_cnt = seq2[42:50].count('A')

                if len1 >= 49 and len2 >= 49 and polyT_cnt>=7 and polyA_cnt>=7:

                    edit_pass1, edit1 = UP_edit_pass(seq1,max_dist)
                    edit_pass2, edit2 = UP_edit_pass(seq2,max_dist)

                    seq_counter(anchor_edits_dict,edit1)
                    seq_counter(target_edits_dict,edit2)

                    if edit_pass1 and edit_pass2:

                        a_bc, a_umi = seq_slice(seq1)
                        t_bc, t_umi = seq_slice(seq2)

                        seq_counter(anchors_umi_len_dict, len(a_umi))
                        seq_counter(targets_umi_len_dict, len(t_umi))

                        quad_dict_store(anchors_dict, a_bc, a_umi)
                        quad_dict_store(targets_dict, t_bc, t_umi)

                    if i>N_read_extract and limit: break
            
            else:
                
                edit_pass1, edit1 = UP_edit_pass(seq1,max_dist)
                edit_pass2, edit2 = UP_edit_pass(seq2,max_dist)

                seq_counter(anchor_edits_dict,edit1)
                seq_counter(target_edits_dict,edit2)

                if edit_pass1 and edit_pass2:

                    a_bc, a_umi = seq_slice(seq1)
                    t_bc, t_umi = seq_slice(seq2)

                    seq_counter(anchors_umi_len_dict, len(a_umi))
                    seq_counter(targets_umi_len_dict, len(t_umi))

                    quad_dict_store(anchors_dict, a_bc, a_umi)
                    quad_dict_store(targets_dict, t_bc, t_umi)

                if i>N_read_extract and limit: break
            
    with open(anchor_edits_json, 'w') as json_file:
        json.dump(anchor_edits_dict, json_file)
    with open(target_edits_json, 'w') as json_file:
        json.dump(target_edits_dict, json_file)
        
    with open(anchors_umi_len_json, 'w') as json_file:
        json.dump(anchors_umi_len_dict, json_file)
    with open(targets_umi_len_json, 'w') as json_file:
        json.dump(targets_umi_len_dict, json_file)
        
    with open(anchors_json, 'w') as json_file:
        json.dump(anchors_dict, json_file)
    with open(targets_json, 'w') as json_file:
        json.dump(targets_dict, json_file)

def extract_quad_dict(indir,sample,part,limit):
    
    i = 0; max_dist = 2; quad_dict = {}
    
    R1_fastq = f'{indir}/{sample}/split/{sample}_R1.part_{part}.fastq'
    R2_fastq = f'{indir}/{sample}/split/{sample}_R2.part_{part}.fastq'
    
    quads_json = f'{indir}/{sample}/split/{sample}.part_{part}_quads.json'
    
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
            
            a_bc, a_umi = seq_slice(seq1)
            t_bc, t_umi = seq_slice(seq2)

            if a_bc in a_dict and t_bc in t_dict:
                edit_pass1, edit1 = UP_edit_pass(seq1,max_dist)
                edit_pass2, edit2 = UP_edit_pass(seq2,max_dist)
                
                if edit_pass1 and edit_pass2:
                    quad_dict_store(quad_dict,a_bc,[a_umi,t_bc])
                    
                    if i>N_read_extract and limit: break
                
    with open(quads_json, 'w') as json_file:
        json.dump(quad_dict, json_file)

def aggregate_stat_dicts(indir,sample,position): 
    
    dir_split=f'{indir}/{sample}/split/'
    files=os.listdir(dir_split)
    jsons = sorted([f for f in files if f'{position}.json' in f])
    
    agg_read_csv=f'{indir}/{sample}/{sample}_agg_cnt_{position}.csv'
    
    if os.path.isfile(agg_read_csv):
        print(agg_read_csv,' exists, skip')
        return
    
    data_agg={}
    for i in tqdm(range(len(jsons))):
        with open(f'{dir_split}{jsons[i]}', 'r') as json_file:
            data_sub = json.load(json_file)
            for k in data_sub:
                if data_agg.get(k) is not None:
                    data_agg[k] += data_sub[k]
                else:
                    data_agg[k]=data_sub[k]
                    
    pd.Series(data_agg).to_csv(agg_read_csv)
    
def aggregate_dicts(indir,sample,position): 
    
    dir_split=f'{indir}/{sample}/split/'
    files=os.listdir(dir_split)
    jsons = sorted([f for f in files if f'{position}.json' in f])
    
    agg_read_csv=f'{indir}/{sample}/{sample}_agg_read_cnt_{position}.csv'
    
    if os.path.isfile(agg_read_csv):
        print(agg_read_csv,' exists, skip')
        return
    
    data_agg = defaultdict(list)
    
    for i in tqdm(range(len(jsons))):
        with open(f'{dir_split}{jsons[i]}', 'r') as json_file:
            data_sub = json.load(json_file)
            print(jsons[i],len(data_sub))
            for key, value in data_sub.items():
                data_agg[key].extend(value)
            
    read_dict = {}
    umi_dict = {}
    total_reads = 0
    for k in tqdm(data_agg):
        reads = len(data_agg[k])
        total_reads += reads
        if reads >= 3:
            umis = len(set(data_agg[k]))
            if umis >= 2:
                read_dict[k] = reads
                umi_dict[k] = umis
            
    print(f'Total Reads Extracted in {position} = {total_reads/1e6}m')
    
    umi_cnt=pd.Series(umi_dict)
    read_cnt=pd.Series(read_dict)
    
    read_cnt.to_csv(agg_read_csv)
    umi_cnt.to_csv(agg_read_csv.replace('read','umi'))
    
def whitelist_rankplot(indir,sample,position,qc_pdfs,max_expected_barcodes=100000):
    
    read_cnt=pd.read_csv(f'{indir}/{sample}/{sample}_agg_read_cnt_{position}.csv')
    umi_cnt=pd.read_csv(f'{indir}/{sample}/{sample}_agg_umi_cnt_{position}.csv')
    umi_cnt.columns=['bc','umi_cnt']
    read_cnt.columns=['bc','read_cnt']
    agg_bcs=pd.merge(umi_cnt,read_cnt,left_on='bc',right_on='bc',how='inner')
    agg_bcs['log10_read_cnt']=np.log10(agg_bcs['read_cnt'])
    agg_bcs['log10_umi_cnt']=np.log10(agg_bcs['umi_cnt'])
    agg_bcs['dup_rate']=agg_bcs['read_cnt']/agg_bcs['umi_cnt']
    agg_bcs=agg_bcs.sort_values(by='umi_cnt',ascending=False)
    
    sub=agg_bcs.iloc[100:max_expected_barcodes].copy()  # select top max_bc except first 100
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

def save_barcode_batch_json(indir,sample):
    
    position='quads'
    
    a_white = pd.read_csv(f'{indir}/{sample}/{sample}_anchors_wl.csv.gz',index_col=1)#['bc']

    dir_split = f'{indir}/{sample}/split/'
    files = os.listdir(dir_split)
    jsons = sorted([f for f in files if f'{position}.json' in f and 'batch_' not in f])
    
    #jsons = jsons#[:subset]
    print(len(jsons))

    data_agg = {}
    
    sub_batch_N = int(len(a_white.index)/30000)+1

    anchors_split = np.array_split(sorted(a_white.index), sub_batch_N)  
    
    for i in tqdm(range(len(jsons))):
        
        part_json_file = f'{dir_split}{jsons[i]}'
        
        batch = str(sub_batch_N).zfill(3)
        batch_json = part_json_file.replace('quads.json',f'batch_{batch}_quads.json')
        if os.path.isfile(batch_json):
            print(batch_json,' exists, skip')
            continue
        
        with open(part_json_file, 'r') as json_file:
            
            data_sub = json.load(json_file)
            print(jsons[i],len(data_sub))

            for j in range(sub_batch_N):
                
                batch = str(j+1).zfill(3)
                batch_json = part_json_file.replace('quads.json',f'batch_{batch}_quads.json')
                
                sub_agg={}
                for a in anchors_split[j]:
                    if a in data_sub:
                        sub_agg[a] = data_sub[a]
              
                print(batch_json,j,len(sub_agg))
                with open(batch_json, 'w') as json_file:
                    json.dump(sub_agg, json_file)

def aggregate_barcode_batches(indir,sample):
    
    dir_split = f'{indir}/{sample}/split/'
    files = os.listdir(dir_split)
    
    jsons = sorted([f for f in files if 'batch_' in f and f'quads.json' in f])
    
    all_batches = np.unique([json.split('batch')[1] for json in jsons])
    
    for b in all_batches:
        
        batch_jsons = ([json for json in jsons if b in json])
        agg_batch_json = f'{sample}.batch{b}'
        agg_batch_json_file = f'{dir_split}{agg_batch_json}'
        if os.path.isfile(agg_batch_json_file):
            print(agg_batch_json_file,' exists, skip')
            continue
        data_agg = {}
        for p in batch_jsons:
            sub_parts_of_batch = f'{dir_split}{p}'
            with open(sub_parts_of_batch, 'r') as json_file:
                data_sub = json.load(json_file)
                print(p,len(data_sub))
                for k in data_sub:
                    if data_agg.get(k) is not None:
                        data_agg[k].extend(data_sub[k])
                    else:
                        data_agg[k] = data_sub[k]
                        
        with open(agg_batch_json_file, 'w') as json_file:
                json.dump(data_agg, json_file)

def make_count_sparse_mtx_batch(indir, sample, batch, threshold=0):
    
    batch = str(batch).zfill(3)
    adata_file = f'{indir}/{sample}/{sample}_counts_b_{batch}.h5ad'
    
    if os.path.isfile(adata_file):
        print(adata_file,' exists, skip')
        return
   
    batch_json = f'{indir}/{sample}/split/{sample}.batch_{batch}_quads.json'

    with open(batch_json, 'r') as json_file:
        data_agg = json.load(json_file)
        
    t_white = pd.read_csv(f'{indir}/{sample}/{sample}_targets_wl.csv.gz')['bc']#,index_col=1)#['bc']
    t_white = t_white.reset_index()
    t_white = t_white.set_index('bc')
    
    a_white = list(data_agg.keys())

    rows_idx = []
    cols_idx = []
    row_col_values = []
    
    for a_idx, a_bc in enumerate(tqdm(a_white)):
        #print(a_idx)
        if a_bc not in data_agg:
            print(f'{a_bc} not in data_agg')
            continue

        umi_tbc = data_agg[a_bc]
        umi_bc_dic = {}
        for a in umi_tbc:
            if len(a[0])==9:
                if umi_bc_dic.get(a[0]) is not None:
                    umi_bc_dic[a[0]].append(a[1])
                else:
                    umi_bc_dic[a[0]] = [a[1]]

        t_bc_cnt={}
        
        for k in umi_bc_dic:
            umi_reads = len(umi_bc_dic[k])
            if umi_reads > threshold:
                uni_t_bc=set(umi_bc_dic[k])
                if len(uni_t_bc) > 1:
                    bcs, cnts = np.unique(umi_bc_dic[k],return_counts=True)
                    if np.max(cnts/umi_reads) > .74: # the concordance to accept the UMI 2/2, 3/3, 3/4, 4/5 or better 
                        t_bc = bcs[np.argmax(cnts/umi_reads)]
                        seq_counter(t_bc_cnt,t_bc)
                else:
                    t_bc = list(uni_t_bc)[0]
                    seq_counter(t_bc_cnt,t_bc)
        
        rows_idx.extend((np.ones(len(t_bc_cnt),dtype=int)*a_idx).tolist())
        cols_idx.extend(t_white.loc[t_bc_cnt.keys()]['index'].tolist())
        row_col_values.extend(list(t_bc_cnt.values()))
        
    csr = csr_matrix((row_col_values, (rows_idx, cols_idx)), shape=(len(a_white), len(t_white)))
    adata = AnnData(csr,dtype='float32')
    adata.var.index = t_white.index
    adata.obs.index = a_white
    
    adata.write_h5ad(adata_file,compression='gzip')
    

def write_fastq_pair_clean(R1_clean, R2_clean, r1, r2, bcs_dict, r1_polyA_cnt_dict, r1_polyT_cnt_dict):
    
    seq1 = r1.sequence
    seq2 = r2.sequence
    
    bc, umi = seq_slice(r1.sequence)
    bc_q, umi_q = seq_slice(r1.quality)
    
    polyT_cnt = seq1[42:50].count('T')
    polyA_cnt = seq1[42:50].count('A')
    seq_counter(r1_polyA_cnt_dict,polyA_cnt)
    seq_counter(r1_polyT_cnt_dict,polyT_cnt)
    """
    polyT_cnt_r2 = seq2.count('T') / len(seq2)
    polyA_cnt_r2 = seq2.count('A') / len(seq2)
    polyG_cnt_r2 = seq2.count('G') / len(seq2)
    
    umi_polyT = umi.count('T')
    
    if umi_polyT<=7 and polyT_cnt_r2<=.6 and polyA_cnt_r2<=.6 and polyG_cnt_r2<=.6:
    """

    quad_dict_store(bcs_dict, bc, umi)

    R1_clean.write(f'@{r1.name}\n')
    R1_clean.write(f'{bc+umi}\n')
    R1_clean.write('+\n')
    R1_clean.write(f'{bc_q+umi_q}\n')

    R2_clean.write(f'@{r2.name}\n')
    R2_clean.write(f'{seq2}\n')
    R2_clean.write('+\n')
    R2_clean.write(f'{r2.quality}\n')
    
def write_fastq_pair_recon(R1_recon, R2_recon, r1, r2):
    
    R1_recon.write(f'@{r1.name}\n')
    R1_recon.write(f'{r1.sequence}\n')
    R1_recon.write('+\n')
    R1_recon.write(f'{r1.quality}\n')

    R2_recon.write(f'@{r2.name}\n')
    R2_recon.write(f'{r2.sequence}\n')
    R2_recon.write('+\n')
    R2_recon.write(f'{r2.quality}\n')
    
def UP_edit_R2(read_seq,max_dist):
    
    edit = edlib.align(UP_seq,read_seq,'HW','distance',max_dist)
    ed_dist = edit['editDistance']
    boolean_fail = ed_dist <0
    
    return(boolean_fail, ed_dist)

def extract_clean_fastq(indir,sample,part,limit):
    
    i = 0
    max_dist = 3
    
    R1_fastq = f'{indir}/{sample}/split/{sample}_R1.part_{part}.fastq'
    R2_fastq = f'{indir}/{sample}/split/{sample}_R2.part_{part}.fastq'
    
    R1_fastq_clean = f'{indir}/{sample}/split/{sample}_R1.part_{part}_clean.fastq'
    R2_fastq_clean = f'{indir}/{sample}/split/{sample}_R2.part_{part}_clean.fastq'
    
    R1_fastq_recon = f'{indir}/{sample}/split/{sample}_R1.part_{part}_recon.fastq'
    R2_fastq_recon = f'{indir}/{sample}/split/{sample}_R2.part_{part}_recon.fastq'
    
    bcs_json = f'{indir}/{sample}/split/{sample}.part_{part}_bcs.json'
    r1_edits_json = f'{indir}/{sample}/split/{sample}.part_{part}_r1_edits.json'
    r1_polyA_cnt_json = f'{indir}/{sample}/split/{sample}.part_{part}_r1_polyA_cnt.json'
    r1_polyT_cnt_json = f'{indir}/{sample}/split/{sample}.part_{part}_r1_polyT_cnt.json'
    
    bcs_dict = {}
    r1_edits_dict = {}
    r1_polyA_cnt_dict = {}
    r1_polyT_cnt_dict = {}

    if os.path.isfile(bcs_json):
        print(bcs_json,' exists, skip')
        return
    
    R1_clean = open(R1_fastq_clean, 'w')
    R2_clean = open(R2_fastq_clean, 'w')
    
    R1_recon = open(R1_fastq_recon, 'w')
    R2_recon = open(R2_fastq_recon, 'w')
    
    with pysam.FastxFile(R1_fastq) as R1, pysam.FastxFile(R2_fastq) as R2:
        for r1, r2 in tqdm(zip(R1, R2)):
            
            i += 1
            
            seq1 = r1.sequence
            seq2 = r2.sequence
            
            len1 = len(seq1)
            len2 = len(seq2)

            if len1 >= 43 and len2 >= 43:
                
                if seq1.count('N')<=5:
                    
                    r1_UP = seq1[8:26] == UP_seq
                    r2_UP = seq2[8:26] == UP_seq

                    if not r2_UP:
                        #edit_fail2, edit2 = UP_edit_R2(seq2,4)
                        #if edit_fail2:
                        if r1_UP:
                            seq_counter(r1_edits_dict,0)
                            write_fastq_pair_clean(R1_clean,R2_clean,r1,r2,bcs_dict,r1_polyA_cnt_dict,r1_polyT_cnt_dict)
                        else:
                            edit_pass1, edit1 = UP_edit_pass(seq1,max_dist)
                            if edit_pass1:
                                seq_counter(r1_edits_dict,edit1)
                                write_fastq_pair_clean(R1_clean,R2_clean,r1,r2,bcs_dict,r1_polyA_cnt_dict,r1_polyT_cnt_dict)
                    else:
                        if r1_UP:
                            write_fastq_pair_recon(R1_recon,R2_recon,r1,r2)

            if i>N_read_extract and limit: break
            
    R1_clean.close()
    R2_clean.close()
    
    R1_recon.close()
    R2_recon.close()
    
    with open(bcs_json, 'w') as json_file:
        json.dump(bcs_dict, json_file)
        
    with open(r1_edits_json, 'w') as json_file:
        json.dump(r1_edits_dict, json_file)
        
    with open(r1_polyA_cnt_json, 'w') as json_file:
        json.dump(r1_polyA_cnt_dict, json_file)
        
    with open(r1_polyT_cnt_json, 'w') as json_file:
        json.dump(r1_polyT_cnt_dict, json_file)
        