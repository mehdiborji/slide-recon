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
parser.add_argument('-ma', '--max_anchors', type=int)
parser.add_argument('-mt', '--max_targets', type=int)

#parser.add_argument('--outdir', type=str)
#parser.add_argument('--barcodes', type=str)
#parser.add_argument('--split', default=False, action='store_true')
#parser.add_argument('--mode', type=str)

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
limit = args.limit
max_anchors = args.max_anchors
max_targets = args.max_targets

#outdir = args.outdir
#barcodes = args.barcodes
#split = args.split
#mode = args.mode

#if __name__ == '__main__':

######################################################

bc_umi_utils.unzip_split_fastq(indir,sample,cores)

#split_fastq(indir,sample,cores)

######################################################

args = bc_umi_utils.find_sub_fastq_pairs(indir,sample,limit)
[print(a) for a in args]

#pool = Pool(int(cores))
pool = Pool(8)

results = pool.starmap(bc_umi_utils.extract_bc_umi_dict, args)
pool.close()
pool.join()

######################################################

bc_umi_utils.aggregate_dicts(indir,sample,'anchors')
bc_umi_utils.aggregate_dicts(indir,sample,'targets')

bc_umi_utils.aggregate_stat_dicts(indir,sample,'anchor_edits')
bc_umi_utils.aggregate_stat_dicts(indir,sample,'target_edits')
bc_umi_utils.aggregate_stat_dicts(indir,sample,'anchors_umi_len')
bc_umi_utils.aggregate_stat_dicts(indir,sample,'targets_umi_len')


######################################################
qc_pdf_file = f'{indir}/{sample}/{sample}_QC.pdf'

if os.path.isfile(qc_pdf_file):
    print(qc_pdf_file,' exists, skip')

else:
    qc_pdfs = PdfPages(qc_pdf_file)
    bc_umi_utils.whitelist_rankplot(indir,sample,'anchors',qc_pdfs,max_anchors)
    bc_umi_utils.whitelist_rankplot(indir,sample,'targets',qc_pdfs,max_targets)
    qc_pdfs.close()
######################################################
args=[(indir, sample, i, limit) for i in range(1, int(cores)+1)]
[print(a) for a in args]

#pool = Pool(int(cores))
pool = Pool(8)

results = pool.starmap(bc_umi_utils.extract_quad_dict, args)
pool.close()
pool.join()
######################################################

#for s in [1, 2, 3, 4, 6, 8, 10, 13, 16]:
for s in [16, 8, 4]:
    bc_umi_utils.make_count_mtx(indir, sample, subset = s, threshold = 0)

######################################################