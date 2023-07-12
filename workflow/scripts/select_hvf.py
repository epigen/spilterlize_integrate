#!/bin/env python

#### select and plot highly variable features (HVF) ####

#### libraries
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statsmodels import robust
import warnings
warnings.filterwarnings("ignore")

# helper functions
def calculate_dispersion(data):
    """
    Calculate the dispersion for each feature across all samples.
    Dispersion is calculated as standard deviation.
    """
    mean = data.mean(axis=1)
    mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    dispersion = data.std(axis=1) # var/mean
    return mean, dispersion

def normalize_dispersion(df, bins_n):
    """
    Normalize the dispersion of each feature by subtracting the median dispersion of its bin and dividing by the MAD of its bin.
    Features are binned based on their means.
    """
    # calcuate the bins
    bins = np.r_[
        -np.inf,
        np.percentile(df['means'], np.arange(0, 100, 100/bins_n)),
        np.inf
    ]
    
    # apply bins
    df['mean_bin'] = pd.cut(df['means'], bins)
    
    # group by means and determine median dispersion
    disp_grouped = df.groupby('mean_bin')['dispersions']
    disp_median_bin = disp_grouped.median()
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # calculate median absolute deviation (mad; robust measure of variance) per group
        disp_mad_bin = disp_grouped.apply(robust.mad)
        # determine "normalized" dispersion via (dispersion-dispersion_group_median)/dispersion_mad
        df['dispersions_norm'] = (df['dispersions'].values
            - disp_median_bin[df['mean_bin'].values].values
            ) / disp_mad_bin[df['mean_bin'].values].values
    return df['dispersions_norm'].values.astype('float32')

    
def visualize_results(df, hvf_idx, bins_n):
    """
    Visualize the results by plotting a histogram of the dispersion values before and after normalization,
    a scatter plot of mean expression levels versus dispersion, and a plot of ranked dispersion.
    The selected features are highlighted in red in the plots.
    """
    plt.figure(figsize=(8,8))

    # Histogram of dispersion values before and after normalization
    plt.subplot(2, 2, 1)
    sns.histplot(df['dispersions'], label='All features')
    sns.histplot(df.loc[hvf_idx, 'dispersions'], color='red', label='Selected HVF')
    plt.xlabel('Dispersion')
    plt.title('Dispersion Before Normalization')
    plt.legend()

    plt.subplot(2, 2, 2)
    sns.histplot(df['dispersions_norm'], label='All features')
    sns.histplot(df.loc[hvf_idx, 'dispersions_norm'], color='red', label='Selected HVF')
    plt.xlabel('Normalized Dispersion')
    plt.title('Dispersion After Normalization')
    plt.legend()

    # Scatter plot of mean expression levels versus normalized dispersion
    plt.subplot(2, 2, 3)
    sns.scatterplot(data=df, x="means", y="dispersions_norm", rasterized=True, label='All features')
    sns.scatterplot(data=df.loc[hvf_idx,:], x="means", y="dispersions_norm", rasterized=True, color='red', label='Selected HVF')
    plt.xlabel('Mean')
    plt.ylabel('Normalized Dispersion')
    plt.title('Mean to Normalized Dispersion')
    plt.legend()

    # Plot of ranked normalized dispersion
    plt.subplot(2, 2, 4)
    df["rank"] = df['dispersions_norm'].rank(ascending=False)
    sns.scatterplot(data=df, x="rank", y="dispersions_norm", rasterized=True, label='All features')
    sns.scatterplot(data=df.loc[hvf_idx,:], x="rank", y="dispersions_norm", rasterized=True, color='red', label='Selected HVF')
    plt.xlabel('Rank')
    plt.ylabel('Normalized Dispersion')
    plt.title('Ranked Dispersion of Features')
    plt.legend()

    plt.gca().set_facecolor('white')
    plt.suptitle('Selection of the top {} highly variable features (HVFs) by binned normalized dispersion'.format(len(hvf_idx)), fontsize=12)
    
    plt.tight_layout()
    plt.savefig(fname=plot_path, format="png", dpi=300, bbox_inches="tight")


#### configurations

# input
data_path = snakemake.input["data"]

# output
result_path = snakemake.output["hvf_data"]
plot_path = snakemake.output["hvf_plot"]

# parameters
hvf_percentage = snakemake.config["hvf_parameters"]["top_percentage"]
bins_n = snakemake.config["hvf_parameters"]["bins_n"]

# load data
data = pd.read_csv(data_path, index_col=0)

# determine number of HVFs
top_n = int((data.shape[0]/100)*hvf_percentage)

# determine normalized dispersion
df = pd.DataFrame()
df['means'], df['dispersions'] = calculate_dispersion(data)
df['dispersions_norm'] = normalize_dispersion(df, bins_n)

# subset & save data
hvf_idx = df.index[df['dispersions_norm'].rank(ascending=False)<=top_n]
hvf_data = data.loc[hvf_idx, :]
hvf_data.to_csv(result_path)

# visualize results
visualize_results(df, hvf_idx, bins_n)
