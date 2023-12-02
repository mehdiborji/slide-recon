"""
plot diffustion stats based on ground truth
input pos and matched reads
output diffusion plot
"""

import os
import random
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.backends.backend_pdf import PdfPages


def plot_cnt_distribution(bead_all, bead_type):
    plt.figure(figsize=(8,6))
    plt.hist(np.log10(bead_all['total_cnt']),100)
    plt.xlabel('log10(total count)')
    plt.ylabel('number of '+bead_type)
    plt.title(bead_type+' total count distribution ({}), median={}'.format(len(bead_all), bead_all['total_cnt'].median()))
    plt.savefig(os.path.join(out_dir,bead_type+'_cnt_distribution.png'),dpi=300)
    plt.close()


def plot_diffusion(pos, match_sum, bead_all, bead_type, base_type):
    match_pos = pd.merge(match_sum, pos, left_on=base_type, right_on='barcode')
    bead_all = bead_all.sort_values(by='total_cnt', ascending=False)
    bi = range(0, 15000, 300)
    pdf = PdfPages(os.path.join(out_dir, bead_type+'_diffusion.pdf'))
    for b_n in bi:
        b = bead_all.iloc[b_n, 0]
        b_match = match_pos.loc[match_pos[bead_type] == b,]
        b_match = b_match.sort_values(by='cnt', ascending=True)
        plt.figure(figsize=(7,5))
        plt.scatter(b_match['xcoord'], b_match['ycoord'], c=np.log(b_match['cnt']+1), alpha=0.8, s=10)
        #plt.scatter(pos.loc[pos['barcode']==b,'xcoord'],pos.loc[pos['barcode']==b,'ycoord'],c='red',marker='.', s=4)
        plt.title(b+' ({})'.format(b_match['cnt'].sum()))
        plt.xlim(0,5000)
        plt.ylim(0,5000)
        cbar = plt.colorbar(pad=0.05, shrink=0.6)
        cbar.set_label('log1p')
        pdf.savefig()
        plt.close()
    pdf.close()


def plot_1d_diffusion(pos, match_sum, bead_all, bead_type, base_type):
    match_pos = pd.merge(match_sum, pos, left_on=base_type, right_on='barcode')
    bi = [random.randint(0, len(bead_all)) for _ in range(100)]
    plt.figure(figsize=(8,6))
    for b_n in bi:
        b = bead_all.iloc[b_n, 0]
        b_match = match_pos.loc[match_pos[bead_type] == b,]
        b_x = np.repeat(b_match['xcoord'],b_match['cnt'])
        b_x -= np.mean(b_x)

        # KDE fitting
        kde = gaussian_kde(b_x)
        x_values = np.linspace(min(b_x), max(b_x), 1000)
        densities = kde.evaluate(x_values)
        plt.plot(x_values, densities, alpha=0.7)

    plt.title(bead_type+' x diffusion')
    plt.xlabel('x')
    plt.ylabel('density')
    plt.savefig(os.path.join(out_dir,bead_type+'_x_diffusion.png'),dpi=300)
    plt.close()


def plot_kde_average(pos, match_sum, bead_all, bead_type, base_type):
    match_pos = pd.merge(match_sum, pos, left_on=base_type, right_on='barcode')
    x_values = np.linspace(-1500, 1500, 1000)
    densities_all = []
    bead_all = bead_all.sort_values(by='total_cnt', ascending=False)
    bi = range(0, 3000)
    for b_n in bi:
        if b_n % 1000 == 0:
            print(b_n)
        b = bead_all.iloc[b_n, 0]
        b_match = match_pos.loc[match_pos[bead_type] == b,]
        if len(b_match) > 4:
            b_x = np.repeat(b_match['xcoord'],b_match['cnt'])
            b_x -= np.mean(b_x)
            kde = gaussian_kde(b_x)
            densities = kde.evaluate(x_values)
            densities_all.append(np.array(densities))
    
    densities_average = np.mean(np.array(densities_all), axis=0)
    plt.figure(figsize=(8,6))
    plt.plot(x_values, densities_average)
    plt.title(bead_type+' x average KDE')
    plt.xlabel('x')
    plt.ylabel('density')
    plt.savefig(os.path.join(out_dir,bead_type+'_x_average_kde.png'),dpi=300)
    plt.close()
    
    return densities_average


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
        "-a", "--anchor",
        help="define anchor bead.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-t", "--target",
        help="define target bead.",
        type=str,
        required=True,
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    sample_folder = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    # sample_folder = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    anchor = args.anchor
    target = args.target

    # make output dir'
    out_dir = os.path.join(sample_folder,'fig')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # read pos and match data
    pos = pd.read_csv(os.path.join(sample_folder,args.sample+'_matched_bead_location.csv.gz'))
    match_raw = pd.read_csv(os.path.join(sample_folder,args.sample+'_raw_matched_reads.csv.gz'))
    match_sum = match_raw.groupby(['P5_bc', 'V9_bc']).size().reset_index(name='cnt')
    match_sum.columns = [target, anchor, 'cnt'] #if tags
    # match_sum.columns = [anchor, target, 'cnt'] #if seq

    # summarise a and t
    a_all = match_sum.groupby(anchor)['cnt'].sum().reset_index(name='total_cnt') 
    t_all = match_sum.groupby(target)['cnt'].sum().reset_index(name='total_cnt')

    # plot a and t cnt distribution
    plot_cnt_distribution(a_all, anchor)
    plot_cnt_distribution(t_all, target)

    # a and t diffusion plot
    plot_diffusion(pos, match_sum, a_all, anchor, target)
    plot_diffusion(pos, match_sum, t_all, target, anchor)

    # 1d diffusion plot
    plot_1d_diffusion(pos, match_sum, a_all, anchor, target)
    plot_1d_diffusion(pos, match_sum, t_all, target, anchor)

    # KDE on average
    a_kde = plot_kde_average(pos, match_sum, a_all, anchor, target)
    a_kde = pd.DataFrame(a_kde, columns=['Value'])
    a_kde.to_csv(os.path.join(out_dir, anchor+'_kde.csv'), index=False)
    t_kde = plot_kde_average(pos, match_sum, t_all, target, anchor)
    t_kde = pd.DataFrame(t_kde, columns=['Value'])
    t_kde.to_csv(os.path.join(out_dir, target+'_kde.csv'), index=False)



