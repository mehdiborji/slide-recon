import argparse
from multiprocessing import Pool
import os
import csv
import gzip
import logging
import numpy as np
import scipy.io
import scipy.sparse
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cores", type=str)
parser.add_argument("-i", "--indir", type=str)

args = parser.parse_args()

cores = args.cores
indir = args.indir


def write_sparse_matrix(dge_txt):
    dge_mtx = dge_txt.replace(".txt", "_matrix.mtx.gz")
    dge_genes = dge_txt.replace(".txt", "_features.tsv.gz")
    dge_barcodes = dge_txt.replace(".txt", "_barcodes.tsv.gz")

    print(dge_mtx)
    # print(dge_genes)
    # print(dge_barcodes)

    print("get cols and rows from dge file")
    with open(dge_txt, "rt") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        cols = next(rdr)[1:]
        rows = [r[0] for r in tqdm(rdr)]

    # print('write barcodes file')
    # with gzip.open(dge_barcodes, "wt") as out:
    #    for bc in cols:
    #        print(bc, file=out)

    print("write features (genes) file")
    with gzip.open(dge_genes, "wt") as out:
        for gene in rows:
            print(gene, file=out)

    data = scipy.sparse.dok_matrix((len(rows), len(cols)), dtype=int)
    with open(dge_txt, "rt") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        _ = next(rdr)
        for i, row in tqdm(enumerate(rdr)):
            for j, val in enumerate(row[1:]):
                if val != "0":
                    data[i, j] = int(val)

    # write mtx file
    with gzip.open(dge_mtx, "wb") as out:
        scipy.io.mmwrite(out, data.tocsr())
    """"""


files = os.listdir(indir)
files = sorted([f"{indir}/{f}" for f in files if "header_added" in f])

pool = Pool(int(cores))
results = pool.map(write_sparse_matrix, files)
pool.close()
pool.join()
