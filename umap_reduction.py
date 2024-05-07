import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
import scanpy as sc
import umap

parser = argparse.ArgumentParser()
parser.add_argument("--cores", type=str)
parser.add_argument("--indir", type=str)
parser.add_argument("--sample", type=str)
parser.add_argument("--adata_name", type=str)
parser.add_argument("--metric", type=str)
parser.add_argument("--n_neighbors", type=int)
parser.add_argument("--min_dist", type=float)
parser.add_argument("--spread", type=float)
# parser.add_argument('--subset', type=int)
# parser.add_argument('--threshold', type=int)

args = parser.parse_args()

cores = args.cores
indir = args.indir
sample = args.sample
adata_name = args.adata_name
metric = args.metric
n_neighbors = args.n_neighbors
min_dist = args.min_dist
spread = args.spread
# subset = args.subset
# threshold = args.threshold


def umap_reduce_batches(
    indir, sample, adata_name, metric, n_neighbors, min_dist, spread
):
    epoch_list = np.unique(
        np.logspace(2.7, 4, num=40, endpoint=True, base=10.0, dtype=int, axis=0)
    ).tolist()

    umap_dir = f"{indir}/{sample}/{adata_name}_{metric}_{n_neighbors}_{min_dist}_1"

    if not os.path.exists(umap_dir):
        os.makedirs(umap_dir)
        print(f"{umap_dir} created")
    else:
        print(f"{umap_dir} already exists")

    epoch_zfilled = f"{str(epoch_list[-1]).zfill(6)}"
    last_umap_csv = f"{umap_dir}/{sample}_e_{epoch_zfilled}.csv"

    if os.path.isfile(last_umap_csv):
        print(last_umap_csv, " exists, skip")
        return

    adata = sc.read(f"{indir}/{sample}/{adata_name}.h5ad")
    reducer = umap.UMAP(
        metric=metric,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        low_memory=False,
        n_components=2,
        verbose=True,
        n_epochs=epoch_list,
        # output_dens = True,
        # local_connectivity = 30,
        learning_rate=1,
    )
    embedding = reducer.fit_transform(adata.X)

    for i, e in enumerate(epoch_list):
        epoch_zfilled = f"{str(e).zfill(6)}"
        epoch_umap = f"{umap_dir}/{sample}_e_{epoch_zfilled}.csv"
        np.savetxt(epoch_umap, reducer.embedding_list_[i], delimiter=",")


def get_umap_limits(indir, sample, subset, threshold):
    umap_added_adata = f"{indir}/{sample}/{sample}_counts_filtered_t_{threshold+1}_s_{subset}_umap_added.h5ad"

    adata = sc.read(umap_added_adata)

    adata_epochs = [col.split("_")[1] for col in adata.obs.columns if "_x" in col]

    ep_cols = [col for col in adata.obs.columns if "ep_" in col]
    min_xy = np.min(adata.obs[ep_cols].values)
    max_xy = np.max(adata.obs[ep_cols].values)
    print(min_xy, max_xy)
    crop_coord = [min_xy, max_xy, min_xy, max_xy]
    return (adata_epochs, crop_coord)


def save_umap_epoch_from_adata(indir, sample, epoch, crop_coord, subset, threshold):
    png_dir = f"{indir}/{sample}/pngs_t_{threshold+1}_s_{subset}"
    if not os.path.exists(png_dir):
        os.makedirs(png_dir)
        # print(f'{png_dir} created')
    else:
        pass
        # print(f'{png_dir} already exists')

    epoch_zfilled = f"{str(epoch).zfill(6)}"
    epoch_png = f"{png_dir}/{sample}_umap_e_{epoch_zfilled}.png"
    if os.path.isfile(epoch_png):
        print(epoch_png, " exists, skip")
        return

    umap_added_adata = f"{indir}/{sample}/{sample}_counts_filtered_t_{threshold+1}_s_{subset}_umap_added.h5ad"

    adata = sc.read(umap_added_adata)

    # adata=sc.read(f'{indir}/{sample}/{sample}_counts_filtered_umap_added.h5ad')
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), gridspec_kw={"wspace": 0.01})
    adata.obsm["spatial"] = adata.obs[[f"ep_{epoch}_x", f"ep_{epoch}_y"]].values
    # sc.pl.umap(adata,crop color='log10_mm10_UMI',s=10,frameon=False)
    sc.pl.spatial(
        adata,
        crop_coord=crop_coord,
        color="log10_gex_UMI",
        spot_size=0.3,
        frameon=False,
        ax=ax,
        show=False,
        cmap="plasma",
    )
    ax.set_title(
        f"epoch_{epoch}", fontsize=13
    )  # , fontweight='bold')#,f'{gene}\n{puck}' fontwieght="medium")
    cbar = ax.collections[0].colorbar
    cbar.remove()
    plt.savefig(epoch_png, bbox_inches="tight")


if __name__ == "__main__":
    
    umap_reduce_batches(
        indir, sample, adata_name, metric, n_neighbors, min_dist, spread
    )

    # epoch_list,crop_coord = get_umap_limits(indir,sample,subset,threshold)

    # args=[(indir,sample,epoch,crop_coord,subset,threshold) for epoch in epoch_list]

    # [print(a) for a in args[:10]]
    # pool = Pool(int(cores))
    # results = pool.starmap(save_umap_epoch_from_adata, args)
    # pool.close()
    # pool.join()
