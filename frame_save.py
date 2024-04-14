import os
import numpy as np
import pandas as pd
import argparse
import subprocess
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument('--cores', type=str)
#parser.add_argument('--indir', type=str)
#parser.add_argument('--sample', type=str)
#parser.add_argument('--adata_name', type=str)

args = parser.parse_args()

cores = args.cores
#indir = args.indir
#sample = args.sample
#adata_name = args.adata_name

def get_umap_limits(indir,sample,subset,threshold):
    
    umap_added_adata = f'{indir}/{sample}/{sample}_counts_filtered_t_{threshold+1}_s_{subset}_umap_added.h5ad'
    
    adata=sc.read(umap_added_adata)
    
    adata_epochs=[col.split('_')[1] for col in adata.obs.columns if '_x' in col]    
    
    ep_cols=[col for col in adata.obs.columns if 'ep_' in col]
    min_xy=np.min(adata.obs[ep_cols].values)
    max_xy=np.max(adata.obs[ep_cols].values)
    print(min_xy,max_xy)
    crop_coord=[min_xy,max_xy,min_xy,max_xy]
    return(adata_epochs,crop_coord)


def save_umap_epoch_from_adata(epoch):
    
    input_dir = '/n/scratch/users/m/meb521/recon'
    
    sample = 'H4_2_rec_merge'
    
    png_dir = f'{input_dir}/{sample}/pngs_withepc_withlimits'
    if not os.path.exists(png_dir):
        try:
            os.makedirs(png_dir)
        except:
            print(f'{png_dir} already created')
    else:
        print(f'{png_dir} already exists')

    epoch_png = f'{png_dir}/umap_e_{epoch}.png'
    if os.path.isfile(epoch_png):
        print(epoch_png,' exists, skip')
        return
    
    # get bcs from a old umap xy merged with barcodes and merged with clusters/colors
    
    umap_run_name = f'{sample}_filtered_.05p_log1p_cosine_35_0.3_nt200_d0.00001_mc200'
    umap = pd.read_csv(f'{input_dir}/{sample}/{umap_run_name}/{sample}_e_010000_umap_bc.csv.gz',index_col=0)
    barcodes = umap.index
    
    
    # get cells and pallete from merged barcodes merged with clusters/colors
    
    cells = pd.read_csv(f'{input_dir}/{sample}/H4_2_first_type.csv',index_col=0)
    pallete = pd.read_csv(f'{input_dir}/{sample}/first_type_palette.csv',index_col=0)
    pal_dict = dict(zip(pallete.type,pallete.palette))
    
    
    # get umap for each epoch of new run
    
    umap_run_name = f'{sample}_filtered_.05p_log1p_cosine_35_0.3_1.0'    
  
    umap = pd.read_csv(f'{input_dir}/{sample}/{umap_run_name}/{sample}_e_{epoch}.csv',names=['x','y'])
    umap.index = barcodes
    
    
    cells_merge = cells.merge(umap, how = 'left', left_index=True, right_index=True)
    plt.rcParams["figure.figsize"] = (10, 10.2)
    sns.scatterplot(data=cells_merge,x='x',y='y',hue='first_type',s=1,palette=pal_dict,legend=None)
    plt.title(f'epoch : {int(epoch)}')
    plt.xlim(-19, 33)  # max limits for x-axis derived from all frames
    plt.ylim(-22, 36)  
    plt.axis('off')
    plt.savefig(epoch_png,bbox_inches='tight');
    plt.close()
    
if __name__ == '__main__':
    
    input_dir = '/n/scratch/users/m/meb521/recon'
    sample = 'H4_2_rec_merge'
    files = sorted(os.listdir(f'{input_dir}/{sample}/H4_2_rec_merge_filtered_.05p_log1p_cosine_35_0.3_1.0'))
    epochs = [f.split('merge_e_')[1].split('.csv')[0] for f in files]
    pool = Pool(int(cores))
    results = pool.map(save_umap_epoch_from_adata, epochs)
    pool.close()
    
    
    
    
    
    
    
    
    
    