import argparse
from multiprocessing import Pool
from matplotlib.backends.backend_pdf import PdfPages
import bc_umi_utils
import os

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

bc_umi_utils.split_fastq_by_lines(indir,sample,16e7)

######################################################

parts = bc_umi_utils.find_sub_fastq_parts(indir,sample)
args = [(indir,sample,part,limit) for part in parts]

pool = Pool(int(cores))
results = pool.starmap(bc_umi_utils.extract_clean_fastq, args)
pool.close()
pool.join()

#bc_umi_utils.aggregate_dicts(indir,sample,'bcs')