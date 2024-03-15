###################
# This plotting script takes output from snakemake.smk piepline in 04_snv_tumor_informed
###################

import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import matplotlib.dates as md
import numpy as np

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['font.family'] = ['Arial']

SMALL_SIZE = 5
MEDIUM_SIZE = 6
BIGGER_SIZE = 7

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def combine_snv_w_background(df, outpath, tumor_sample, every_x=5):
    fig, ax1 = plt.subplots(figsize=( 7.20472/2, 7.20472/4 ))
    # Plot the data
    ax1.set_zorder(1)
    ax1.patch.set_visible(False)
    for sample in df['sample'].unique():
        if sample != tumor_sample:
            grouped = df[df['sample'] == sample]
            color = 'grey'
            ax1.plot(grouped['time_bin'], grouped['ratio'], linewidth=0.5, marker='o', markersize = 1.5, color=color)
    tumor = df[df['sample'] == tumor_sample]
    color = 'darkred'
    ax1.plot(tumor['time_bin'], tumor['ratio'], linewidth=0.5, marker='o', markersize = 1.5, color=color)

    locs, _ = plt.xticks()  # Get current x-tick locations
    ## If not fixed, ValueError: The number of FixedLocator locations (13), usually from a call to set_ticks, does not match the number of ticklabels (12).
    minutes = df['minutes'].unique()
    plt.xticks(locs[::every_x], labels=[int(x) for x in minutes[::every_x]], rotation=45)  # Set x-tick locations to every 10th tick

    step = 0.05
    ax1.set_yticks(ticks=np.arange(0, 0.5 + step , step ))
    ax1.set_ylim(0, 0.5)
    # ax1.set_yticks(ticks=np.arange(-0.15, 1 + step , step ), labels = [""]*3 + list(np.arange(-0.15, 1 + step , step ))[3:])
    ax1.set_xlabel('Minute sequencing')
    ax1.set_ylabel('ALT Ratio')
    # Plot the background filled SNV count
    ax2 = ax1.twinx()
    ax2.set_ylabel('SNV count')
    ax2.fill_between(tumor['time_bin'], tumor['REF_sum']+tumor['ALT_sum'], interpolate=True, facecolors="#d8d8d8", edgecolor= "#d8d8d8")
    ax2.fill_between(tumor['time_bin'], tumor['ALT_sum'], interpolate=True, facecolors="#c09696", edgecolor= "#c09696")
    if tumor['sum_snv'].max() < 70:
        step2 = 5
    elif tumor['sum_snv'].max() < 200:
        step2 = 20
    elif  tumor['sum_snv'].max() < 1000:
        step2 = 50
    else:
        step2 = 100
    ax2.set_yticks(ticks=np.arange(0, tumor['sum_snv'].max() + step2 , step2 ))
    ax2.set_ylim(0, tumor['sum_snv'].max() + step2/2 )
    plt.title(tumor_sample)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.savefig(f"{outpath}.png", dpi=300)
    plt.show()

    # plt.show()



if __name__ == '__main__':
    # Define output path and output name
    samples = ["OVCA01", "GCT01","GCT02","OES01","OES02",]

    total_mins = 1200
    max_plotting_minutes = 1200

    for tumor_sample in samples:
        outpath = f"/path/to/NanoRCS/output/output_figures/Fig2D/{tumor_sample}_realtime_combined_{max_plotting_minutes}.pdf"
        path_prefix = "/path/to/NanoRCS/output/processed_data/04_snv_tumor_informed/overlap/40/realtime/"
        TUMOR_df = pd.read_pickle(f"{path_prefix}/{tumor_sample}_{tumor_sample}_SNV_realtime_{total_mins}.pickle.gz")
        HC01_df =  pd.read_pickle(f"{path_prefix}/HC01_{tumor_sample}_SNV_realtime_{total_mins}.pickle.gz")
        HC02_df =  pd.read_pickle(f"{path_prefix}/HC02_{tumor_sample}_SNV_realtime_{total_mins}.pickle.gz")
        HC03_df =  pd.read_pickle(f"{path_prefix}/HC03_{tumor_sample}_SNV_realtime_{total_mins}.pickle.gz")
        # This would fill the remaining bins with previous value if value is empty (For sequencing which is still on-going, but no ALT nor REF found.
        TUMOR_df = TUMOR_df.fillna(method='ffill')
        HC01_df = HC01_df.fillna(method='ffill')
        HC02_df = HC02_df.fillna(method='ffill')
        HC03_df = HC03_df.fillna(method='ffill')

        ## input params
        # set max mins = 60
        # in second ,fixed
        bin_size = 60
        # taken from the realtime.py
        # Every Y minutes a dot
        Y_min = 10
        # Ever X minutes a tick
        X_min = 60
        every_x = X_min // Y_min
        # main
        df123 = pd.concat([TUMOR_df, HC01_df, HC02_df, HC03_df ])
        # update df123
        time_bins = [int(v.rsplit('-',1)[1]) for v in list(df123.index)]
        minutes = np.array(time_bins) / 60
        df123['minutes'] = minutes
        df123 = df123.reset_index()
        # Select the dots to plot
        df_subset_to_plot = df123[df123['minutes'].isin(list( range(0, max_plotting_minutes+Y_min,Y_min )))]
        combine_snv_w_background(df_subset_to_plot, outpath, tumor_sample, every_x = every_x)
