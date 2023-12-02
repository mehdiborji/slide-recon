"""
analysis of fiducial diffusion sequencing result without bead barcode matching
input fastq file
output collasped barcode information
"""

import os
import time
import gzip
import argparse
import numpy as np
import mappy as mp
import editdistance
import matplotlib.pyplot as plt
from collections import Counter
from multiprocessing import Pool
from umi_tools import UMIClusterer
from bead_matching import barcode_matching


def barcode_extract(fq1_file, fq2_file):
    """
    input: fastq file
    output: dict of barcode without matching
    """
    aln_dict = {}
    P5_bc_list = []
    V9_bc_list = []
    alignment_stat = Counter()
    for fq1, fq2 in zip(mp.fastx_read(fq1_file, read_comment=False), mp.fastx_read(fq2_file, read_comment=False)):
        alignment_stat["total_reads"] += 1

        if alignment_stat["total_reads"] % 1000000 == 0:
            print(alignment_stat["total_reads"])
    
        # check read length
        if len(fq1[1])<41 or len(fq2[1])<41:  
            alignment_stat["Read_too_short"] += 1
            continue
        P5_bc = fq1[1][:8]+fq1[1][26:32]
        P5_bumi = fq1[1][33:41] # if P5: 32:40, if V8 or V10: 33:41
        P5_UP = fq1[1][8:26]
        V9_bc = fq2[1][:8]+fq2[1][26:32]
        V9_bumi = fq2[1][33:41]
        V9_UP = fq2[1][8:26]

        # check UP site with mismatch < 3bp
        if editdistance.eval(P5_UP,'TCTTCAGCGTTCCCGAGA')>3 or editdistance.eval(V9_UP,'TCTTCAGCGTTCCCGAGA')>3:
            alignment_stat["UP_not_matched"] += 1
            continue 

        aln_dict.setdefault(P5_bc,[]).append((V9_bc,P5_bumi,V9_bumi)) 
        P5_bc_list.append(P5_bc)
        V9_bc_list.append(V9_bc)
    return aln_dict, alignment_stat, P5_bc_list, V9_bc_list


def umi_collapsing(cnt_dict, max_dist=1):
    """
    input: dict of barcode without matching
    output: list of barcode after collapsing
    """
    start_time = time.time()
    clusterer = UMIClusterer(cluster_method="directional")
    clustered_bc = clusterer(cnt_dict, threshold=max_dist)
    clustering_time = time.time()
    cluster_bc = [bc_group[0].decode('utf-8') for bc_group in clustered_bc]
    end_time = time.time()
    print("Clustering time: {}s".format(clustering_time-start_time))
    print("Dict creation time is: {}s".format(end_time-clustering_time))
    print("Total time is: {}s".format(end_time-start_time))
    return cluster_bc


def bc_collapsing(aln_dict, P5_bc_list, V9_bc_list, min_reads_P5, min_reads_V9, alignment_stat):
    """ 
    input: dict of barcode without matching
    output: dict of barcode after filtering and collapsing
    """
    # filter for reads and collapse to whitelist
    P5_list = [s.encode('utf-8') for s in P5_bc_list]
    P5_dict = dict(Counter(P5_list))
    P5_dict_top = {k: v for k, v in P5_dict.items() if v > min_reads_P5}
    P5_whitelist = umi_collapsing(P5_dict_top)
    print("P5 total {}, after filter {}, whitelist {}".format(len(P5_dict),len(P5_dict_top),len(P5_whitelist)))
    print("read percentage: {}".format(np.sum(list(P5_dict_top.values()))/np.sum(list(P5_dict.values()))))
    V9_list = [s.encode('utf-8') for s in V9_bc_list]
    V9_dict = dict(Counter(V9_list))
    V9_dict_top = {k: v for k, v in V9_dict.items() if v > min_reads_V9}
    V9_whitelist = umi_collapsing(V9_dict_top)
    print("V9 total {}, after filter {}, whitelist {}".format(len(V9_dict),len(V9_dict_top),len(V9_whitelist)))
    print("read percentage: {}".format(np.sum(list(V9_dict_top.values()))/np.sum(list(V9_dict.values()))))

    # match to whitelist
    P5_bc_matching_dict,_,_ = barcode_matching(Counter(P5_whitelist), list(set(P5_bc_list)), max_dist=1)
    V9_bc_matching_dict,_,_ = barcode_matching(Counter(V9_whitelist), list(set(V9_bc_list)), max_dist=1)

    # generate dict with matched bc
    aln_dict_new = {}
    for bc_P5 in aln_dict:
        if bc_P5 in P5_bc_matching_dict:
            for V9 in range(len(aln_dict[bc_P5])):
                bc_V9 = aln_dict[bc_P5][V9][0]
                if bc_V9 in V9_bc_matching_dict:
                    alignment_stat["matched_reads"] += 1
                    aln_dict_new.setdefault(P5_bc_matching_dict[bc_P5],[]).append(
                        (V9_bc_matching_dict[bc_V9],aln_dict[bc_P5][V9][1],aln_dict[bc_P5][V9][2])) 
    print(len(aln_dict_new))
                    
    return aln_dict_new, alignment_stat


def write_blind(aln_dict_new, alignment_stat, sample, out_dir):
    # collapse for reads
    for bc in aln_dict_new:
        tmp = Counter(aln_dict_new[bc])
        aln_dict_new[bc] = tmp

    # write result to csv
    raw_f = gzip.open(os.path.join(out_dir,(sample+"_blind_raw_reads_filtered.csv.gz")),"wb")
    raw_f.write(b'P5_bc,V9_bc,P5_bumi,V9_bumi,reads\n')
    for bc_P5 in aln_dict_new:
        raw_f.write(bytes('\n'.join(['{},{},{},{},{}'.format(
            bc_P5, it[0], it[1], it[2], aln_dict_new[bc_P5][it]) for it in aln_dict_new[bc_P5]])+'\n',"UTF-8"))
    raw_f.close()
    print("Write matched data to {}".format("blind_raw_reads_filtered.csv.gz"))

    with open(os.path.join(out_dir,(sample+"_blind_statistics_filtered.csv")),"w") as f:
        f.write("alignment_status,counts\n")
        for aln_stat in alignment_stat:
            f.write("{},{}\n".format(aln_stat, alignment_stat[aln_stat]) )


def get_args():
    parser = argparse.ArgumentParser(description='Process recon seq data.')
    parser.add_argument("-d", "--date",
        help="input experiment data.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s", "--sample",
        help="input sample id.",
        type=str,
        required=True,
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()

    fq_dir = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    # fq_dir = os.path.join("/Volumes/broad_thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    # fq_dir = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    fq_files = os.listdir(fq_dir)
    fq1_file = os.path.join(fq_dir,[it for it in fq_files if args.sample in it and "R1" in it][0])
    fq2_file = os.path.join(fq_dir,[it for it in fq_files if args.sample in it and "R2" in it][0])

    # make output dir'
    out_dir = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample+'_30_100')
    # out_dir = os.path.join("/Volumes/broad_thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    # out_dir = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # bacode blind
    aln_dict, stat, P5_bc_list, V9_bc_list = barcode_extract(fq1_file, fq2_file)
    print("barcode extracted")

    # plot bc rank
    P5_bc_dict = Counter(P5_bc_list).most_common()
    P5_bc_rank = [cnt for _, cnt in P5_bc_dict]
    x_val = range(1,len(P5_bc_rank)+1)
    plt.plot(x_val, P5_bc_rank, marker='o', markersize=2, linestyle='None')
    plt.xscale('log') 
    plt.yscale('log')
    plt.xlabel('bc Rank')
    plt.ylabel('read Counts')
    plt.title('Barcode Rank Plot')
    plt.savefig(os.path.join(out_dir,'P5_barcode_rank_plot_loglog.png'),dpi=300)
    plt.close()

    V9_bc_dict = Counter(V9_bc_list).most_common()
    V9_bc_rank = [cnt for _, cnt in V9_bc_dict]
    x_val = range(1,len(V9_bc_rank)+1)
    plt.plot(x_val, V9_bc_rank, marker='o', markersize=2, linestyle='None')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('bc Rank')
    plt.ylabel('read Counts')
    plt.title('Barcode Rank Plot')
    plt.savefig(os.path.join(out_dir,'V9_barcode_rank_plot_loglog.png'),dpi=300)
    plt.close()

    aln_dict_new, stat_new = bc_collapsing(aln_dict, P5_bc_list, V9_bc_list, min_reads_P5=30, min_reads_V9=100, alignment_stat = stat)
    write_blind(aln_dict_new, stat_new, args.sample, out_dir)
