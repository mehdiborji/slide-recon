"""
analysis of fiducial diffusion sequencing result
input fastq file
output matching information and also barcode ground truth
"""

import os
import gzip
import argparse
import mappy as mp
import editdistance
from collections import Counter
from scipy.spatial.distance import hamming
from bead_matching import get_barcode_position, barcode_matching

def match_groundtruth(fq1_file, fq2_file, bc_pos_dict):
    """
    input: fastq file
    output: dict of matched barcode
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
        if len(fq1[1])<40 or len(fq2[1])<41:  
            alignment_stat["Read_too_short"] += 1
            continue
        a_bc = fq1[1][:8]+fq1[1][26:32]
        a_bumi = fq1[1][33:41] # if P5: 32:40, if V8 or V10: 33:41
        a_UP = fq1[1][8:26]
        t_bc = fq2[1][:8]+fq2[1][26:32]
        t_bumi = fq2[1][33:41]
        t_UP = fq2[1][8:26]

        # check UP site with mismatch < 3bp
        if editdistance.eval(list(a_UP),list('TCTTCAGCGTTCCCGAGA'))>3 or editdistance.eval(list(t_UP),list('TCTTCAGCGTTCCCGAGA'))>3:
            alignment_stat["UP_not_matched"] += 1
            continue 

        aln_dict.setdefault(a_bc,[]).append((t_bc,a_bumi,t_bumi)) 
        P5_bc_list.append(a_bc)
        V9_bc_list.append(t_bc) 
    print('read fastq done')

    # match bead barcode
    P5_bc_matching_dict,_,_ = barcode_matching(bc_pos_dict, P5_bc_list, max_dist=1)    
    V9_bc_matching_dict,_,_ = barcode_matching(bc_pos_dict, V9_bc_list, max_dist=1)
    print('barcode matching done')  

    # generate dict with matched bc
    aln_dict_new = {}
    for bc_a in aln_dict:
        if bc_a in P5_bc_matching_dict:
            for t in range(len(aln_dict[bc_a])):
                bc_t = aln_dict[bc_a][t][0]
                if bc_t in V9_bc_matching_dict:
                    alignment_stat["both_bc_matched"] += 1
                    aln_dict_new.setdefault(P5_bc_matching_dict[bc_a],[]).append(
                        (V9_bc_matching_dict[bc_t],aln_dict[bc_a][t][1],aln_dict[bc_a][t][2])) 
    
    return aln_dict_new, alignment_stat


def write_groundtruth(aln_dict_new, alignment_stat, bc_pos_dict, sample, out_dir):
    # write bead location
    loc_file = os.path.join(out_dir,(sample+'_matched_bead_location.csv.gz'))
    with gzip.open(loc_file,"wb") as f:
        f.write(b'barcode,xcoord,ycoord\n')
        for bc in bc_pos_dict:
            f.write(bytes("{},{},{}\n".format(bc, bc_pos_dict[bc][0],bc_pos_dict[bc][1] ),"UTF-8"))

    # collapse for reads
    for bc in aln_dict_new:
        tmp = Counter(aln_dict_new[bc])
        aln_dict_new[bc] = tmp

    # write result to csv
    raw_f = gzip.open(os.path.join(out_dir,(sample+"_raw_matched_reads.csv.gz")),"wb")
    raw_f.write(b'P5_bc,V9_bc,P5_bumi,V9_bumi,reads\n')
    for bc_P5 in aln_dict_new:
        raw_f.write(bytes('\n'.join(['{},{},{},{},{}'.format(
            bc_P5, it[0], it[1], it[2], aln_dict_new[bc_P5][it]) for it in aln_dict_new[bc_P5]])+'\n',"UTF-8"))
    raw_f.close()
    print("Write matched data to {}".format("raw_matched_reads.csv.gz"))

    with open(os.path.join(out_dir,(sample+"_aligment_statistics.csv")),"w") as f:
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
    parser.add_argument(
        "-p", "--puckid",
        help="the puck_id such as Puck_211004_34. will be looking for folder with the same name in /broad/macosko/data/Slideseq/Barcode.",
        type=str,
        required=True,
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()

    puck_dir = os.path.join("/broad/macosko/data/Slideseq/Barcodes",args.puckid)
    # puck_dir = os.path.join("/mnt/macoskolab/data/Slideseq/Barcodes",args.puckid)
    bead_bc_file = os.path.join(puck_dir,"BeadBarcodes.txt")
    bead_pos_file = os.path.join(puck_dir,"BeadLocations.txt")
    bc_pos_dict, bc_list = get_barcode_position(bead_bc_file,bead_pos_file)   

    fq_dir = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    # fq_dir = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    fq_files = os.listdir(fq_dir)
    fq1_file = os.path.join(fq_dir,[it for it in fq_files if args.sample in it and "R1" in it][0])
    fq2_file = os.path.join(fq_dir,[it for it in fq_files if args.sample in it and "R2" in it][0])

    # make output dir'
    out_dir = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    # out_dir = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # matching
    aln_dict, stat = match_groundtruth(fq1_file,fq2_file,bc_pos_dict)
    print("match done")
    write_groundtruth(aln_dict,stat,bc_pos_dict,args.sample,out_dir)



    