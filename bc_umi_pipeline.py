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

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
limit = args.limit
max_anchors = args.max_anchors
max_targets = args.max_targets

######################################################

bc_umi_utils.split_fastq_by_lines(indir,sample,4e6)

######################################################

parts = bc_umi_utils.find_sub_fastq_parts(indir,sample)
args = [(indir,sample,part,limit) for part in parts]

pool = Pool(int(cores))
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

pool = Pool(int(cores))
results = pool.starmap(bc_umi_utils.extract_quad_dict, args)
pool.close()
pool.join()
######################################################

bc_umi_utils.save_barcode_batch_json(indir,sample)
bc_umi_utils.aggregate_barcode_batches(indir,sample)
                                     
######################################################

batches = sorted([f.split('_')[-2] for f in os.listdir(f'{indir}/{sample}/split/') if 'batch' in f and 'part' not in f])
len_batches = len(batches)
args = [(indir, sample, i) for i in range(1, len_batches+1)]
[print(a) for a in args]

pool = Pool(int(cores))
results = pool.starmap(bc_umi_utils.make_count_sparse_mtx_batch, args)
pool.close()
pool.join()