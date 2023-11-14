"""
Reconstuction from h5ad matrix to 2D embedding, compare with ground truth
input h5ad file with anchor by target matrix
output 2D embedding and comparison
"""

import os
import time
import umap
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import editdistance
import colorsys as cs
from collections import Counter
import matplotlib.pyplot as plt
from scipy.spatial import procrustes
from scipy.spatial.distance import cdist

# generate coordinates correlated hsl color
def get_col(coords):
    plot_df = coords.copy()
    my_h = plot_df['xcoord']
    my_h = (my_h - np.min(my_h))
    my_h = 0.2+my_h * (0.8 / np.max(my_h))
    my_s = np.ones(len(my_h)) * 0.5
    my_l = plot_df['ycoord']
    my_l = (my_l - np.min(my_l))
    my_l = 0.2 + my_l * (0.7 / np.max(my_l))
    hls_color = np.column_stack((my_h, my_l, my_s))
    rgb_color = [cs.hls_to_rgb(p[0], p[1], p[2]) for p in hls_color]
    return rgb_color

# get ground truth 
def get_truth(wlist,pos_df):
    w_truth = pos_df[pos_df['barcode'].isin(wlist)]
    w_truth = w_truth.sort_values(by=['barcode'])
    w_col = get_col(w_truth)
    return w_truth, w_col

# procrustes analysis for alignment
def get_pa(pos_truth,pos_recon):
    mtx1, mtx2, disparity = procrustes(pos_truth[['xcoord','ycoord']], pos_recon[['xcoord','ycoord']])
    
    pos_truth_translate = pos_truth.copy()
    pos_truth_translate['xcoord'] = pos_truth_translate['xcoord'] - np.mean(pos_truth_translate['xcoord'])
    pos_truth_translate['ycoord'] = pos_truth_translate['ycoord'] - np.mean(pos_truth_translate['ycoord'])
    
    scaling = np.sqrt(np.trace(np.dot(pos_truth_translate[['xcoord','ycoord']], pos_truth_translate[['xcoord','ycoord']].T)))
    scaling_2 = np.sqrt(np.trace(np.dot(pos_truth_translate[['xcoord','ycoord']], pos_truth_translate[['xcoord','ycoord']].T)))/np.sqrt(np.trace(np.dot(mtx2, mtx2.T)))

    mtx1_scaled = mtx1*scaling*0.72
    mtx2_scaled = mtx2*scaling_2*0.72
    
    dist = pos_truth_translate.copy()
    dist['x'] = mtx1_scaled[:,0] - mtx2_scaled[:,0]
    dist['y'] = mtx1_scaled[:,1] - mtx2_scaled[:,1]
    dist['r'] = np.sqrt(dist['x']**2 + dist['y']**2)
    
    return disparity, mtx1_scaled, mtx2_scaled, dist


# sample info
date = '230927'
sample = 'c58_4'
anchor = 'V9A30'
target = 'V10'

sample_filter_folder = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",date,'filtered_data')
sample_folder = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",date,sample)
folder_files = os.listdir(sample_filter_folder)
adata_file = os.path.join(sample_filter_folder,[it for it in folder_files if sample in it and "h5ad" in it][0])


# filter for count>0 and generate count matrix
adata=sc.read_h5ad(adata_file)
adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.var['n_counts'] = np.array(adata.X.sum(axis=0)).flatten()
adata_filtered = adata[adata.obs['n_counts'] > 0, adata.var['n_counts'] > 0]
counts = adata_filtered.X.toarray()


# make output folder for recon results
out_dir = os.path.join(sample_folder,'recon')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)


# recon by UMAP
time0 = time.time()
reducer = umap.UMAP(metric='cosine',
                    n_neighbors=25, 
                    min_dist=0.99, 
                    low_memory=False, 
                    n_components=2, 
                    # random_state=0, 
                    verbose=True, 
                    n_epochs=2000,
                    # output_dens = True,
                    # local_connectivity = 30,
                    learning_rate = 1)
embedding = reducer.fit_transform(counts)
print('umap time: ', time.time()-time0)

plt.figure(figsize=(6,6))
plt.scatter(embedding[:,0],embedding[:,1],s=1)
plt.title(anchor+' umap ({})'.format(len(embedding)))
# plt.savefig(os.path.join(out_dir,'{}_UMAP.png'.format(anchor)),dpi=300)
# plt.close()
plt.show()


# match to ground truth
pos = pd.read_csv(os.path.join(sample_folder, sample+'_matched_bead_location.csv.gz'))
bc_pos_dict = Counter(pos['barcode'])
a_bc_matching_dict,_,_ = barcode_matching(bc_pos_dict, adata_filtered.obs.index.values, max_dist=1) 
a_truth, a_col = get_truth(list(a_bc_matching_dict.values()), pos)

recon_umap = pd.DataFrame(embedding)
a_umap  = recon_umap.copy()
a_umap[anchor] = adata_filtered.obs.index.values
a_umap.columns = ['xcoord', 'ycoord', anchor]
a_umap_matched = a_umap[a_umap[anchor].isin(a_bc_matching_dict.keys())]
a_umap_matched = a_umap_matched.copy()
a_umap_matched.loc[:,'barcode'] = a_umap_matched[anchor].map(a_bc_matching_dict)
a_umap_matched = a_umap_matched[['xcoord','ycoord','barcode']]
a_umap_matched = a_umap_matched.groupby('barcode').mean().reset_index()
a_umap_matched = a_umap_matched.sort_values(by=['barcode'])

# 2D embedding compared to ground truth
fig, axes = plt.subplots(ncols=2, figsize=(12, 5))
axes[0].scatter(a_truth['xcoord'], a_truth['ycoord'], s=3, c=a_col)
axes[0].set_title('ground truth ({})'.format(len(a_truth)))

axes[1].scatter(a_umap_matched['xcoord'], a_umap_matched['ycoord'], s=3, c=a_col)
axes[1].set_title(anchor+'UMAP ({})'.format(len(a_umap_matched)))
    
plt.show()


# alignment analysis
da, a_truth_pa, a_recon_pa, a_comp = get_pa(a_truth,a_umap_matched)
print('difference: {}'.format(da))

print('mean distances: {}'.format(np.mean(a_comp['r'])))
print('median distances: {}'.format(np.median(a_comp['r'])))

plt.figure(figsize=(6,3))
plt.hist(a_comp['r'],bins=1000)
plt.xlim(0,600)
plt.title(anchor+' umap differece, median={}'.format(np.median(a_comp['r'])))
# plt.savefig(os.path.join(out_dir,'{}_UMAP_PA.png'.format(anchor)),dpi=300)
# plt.close()
plt.show()


# compare of relative error
a_comp.reset_index(drop=True, inplace=True)
indices = a_comp[a_comp['r'] < 160].index.tolist()

a_truth_pa_pairwise = cdist(a_truth_pa[indices,:], a_truth_pa[indices,:])
a_recon_pa_pairwise = cdist(a_recon_pa[indices,:], a_recon_pa[indices,:])
diff = a_recon_pa_pairwise - a_truth_pa_pairwise

intervals = [(i*30+0.001, (i+1)*30) for i in range(100)]

rms_values = []
a_truth_pairwise = a_truth_pa_pairwise.flatten()
diff = diff.flatten()

for interval in intervals:
    mask = (a_truth_pairwise >= interval[0]) & (a_truth_pairwise < interval[1])
    corresponding_values = diff[mask]
    rms = np.sqrt(np.mean(corresponding_values**2))
    rms_values.append(rms)

plt.figure(figsize=(6,4))
plt.plot([i*30 for i in range(100)],rms_values)
plt.xlabel('distance in truth/um')
plt.ylabel('RMS of distance difference/um')
plt.title(anchor+' error at distance (UMAP)')
plt.show()

intervals = [(i*5+0.001, (i+1)*5) for i in np.arange(2,100)]

rms_values = []
a_truth_pairwise = a_truth_pa_pairwise.flatten()
diff = diff.flatten()

for interval in intervals:
    mask = (a_truth_pairwise >= interval[0]) & (a_truth_pairwise < interval[1])
    corresponding_values = diff[mask]
    rms = np.sqrt(np.mean(corresponding_values**2))
    rms_values.append(rms)

print('100 um: ', rms_values[18]/100)
print('200 um: ', rms_values[38]/200)
print('300 um: ', rms_values[58]/300)
print('400 um: ', rms_values[78]/400)

plt.figure(figsize=(6,4))
plt.plot([i*5 for i in np.arange(2,100)],[(rms_values[i-2])/(i*5) for i in np.arange(2,100)], label='recon')
plt.plot([i*5 for i in np.arange(2,100)],[10/(i*5) for i in np.arange(2,100)], label='Slide-seq')
plt.xlabel('distance in truth/um')
plt.legend()
plt.ylabel('RMS of distance difference/distance')
plt.title(anchor+'relative error at distance (UMAP)')
plt.xscale('log') 
plt.show()