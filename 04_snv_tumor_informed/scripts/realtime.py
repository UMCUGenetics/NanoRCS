import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import matplotlib.dates as md
import numpy as np
import os
import argparse





def get_read_id(bam):
    read_id_list = []
    with pysam.AlignmentFile(bam, 'r') as f:
        for read in f:
            read_id = read.query_name.split('_')[0]
            read_id_list.append(read_id)
    if len(read_id_list) == 0:
        print(f"Warning: This bam file has no reads. f{bam}")

    return read_id_list


def get_timestamp_of_reads(df, sample_name, snv_type, id_list=None):
    if id_list != None:
        df = df[df['read_id'].isin(id_list)].copy()
    df['time'] = df['start_time'] + df['duration']
    df = df[['time', 'read_id']].copy()
    df['sample'] = sample_name
    df['snv_type'] = snv_type
    df['timedelta'] = pd.to_timedelta(df['time'], unit='s')
    return df


def get_read_subset_with_timestamp(summary_df, bam, sample_name, snv_type):
    id_list = get_read_id(bam)
    df = get_timestamp_of_reads(summary_df, sample_name, snv_type, id_list)
    return df


### Choose same bin size as grouped_snv
def plot_cum(df, outpath):
    sample = df['sample'].values[0]
    palette = {
        'ALT': 'red',
        'REF': 'grey'
    }
    # Using Seaborn's ECDF plot
    sns.ecdfplot(data=df, x='timedelta', hue='snv_type', stat="count", palette=palette)

    # Formatting the x-axis ticks to show hours, minutes, seconds
    def format_timedelta(tick_val, _):
        # Convert from nanoseconds to seconds
        seconds = tick_val / 1e9
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"

    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_timedelta))
    plt.xticks(rotation=90)
    plt.xlabel('Relative Time (HH:MM:SS)')
    plt.ylabel('ECDF')
    plt.title(f'{sample} SNV count over time')
    plt.tight_layout()
    plt.savefig(outpath)
    plt.savefig(f"{outpath}.png", dpi =300)
    # plt.show()


def plot_alt_ratio(grouped, bins, outpath):
    sample = grouped['sample'].values[0]
    minutes = np.array(list(bins)[1:]) / 60
    every_x = 5
    #     plt.xticks(locs[::every_x], labels=minutes[::every_x])  # Set x-tick locations to every 10th tick

    fig, ax1 = plt.subplots(figsize=(8, 4))
    # Plot the data
    ax1.set_zorder(1)
    ax1.patch.set_visible(False)
    ax1.plot(grouped['ratio'], marker='o', color='darkred')
    locs, _ = plt.xticks()  # Get current x-tick locations
    ## If not fixed, ValueError: The number of FixedLocator locations (13), usually from a call to set_ticks, does not match the number of ticklabels (12).
    plt.xticks(locs[::every_x], labels=minutes[::every_x])  # Set x-tick locations to every 10th tick
    step = 0.05
    plt.yticks(ticks=np.arange(0, 1 + step , step ))
    ax1.set_xlabel('Minute sequencing')
    ax1.set_ylabel('ALT Ratio')
    # ax1.title('ALT allele ratio over first hour of sequencing')
    # ax1.xticks(rotation=90)  # Optional: rotate x-tick labels for better readability
    # ax1.tick_params(axis='y', labelcolor='b')

    # # Create the secondary y-axis and plot on it
    # ax2 = ax1.twinx()
    # ax2.fill_between(grouped.index, grouped['sum_snv'], interpolate=True, color='orange', alpha=0.3)
    # ##ax2.plot(grouped['sum_snv'], marker = '.', color='orange')
    # ax2.set_ylabel('Total SNV sequenced', color='orange')
    # ## ax2.tick_params(axis='y', labelcolor='r')

    ax2 = ax1.twinx()
    ax2.fill_between(grouped.index, grouped['sum_snv'], interpolate=True, color='grey', alpha=0.3)
    ax2.fill_between(grouped.index, grouped['sum_snv'], interpolate=True, color='grey', alpha=0.3)
    ax2.set_ylabel('Total read', color='darkgreen')

    plt.title(f'{sample} ALT allele ratio over first hour of sequencing')
    plt.xticks(rotation=90)  # Optional: rotate x-tick labels for better readability
    plt.tight_layout()  # Adjust layout for better visualization
    plt.savefig(f"{outpath}/{sample}_ALT_ratio_overtime_throughput.pdf", dpi=300)
    plt.savefig(f"{outpath}/{sample}_ALT_ratio_overtime_throughput.png", dpi=300)
    # plt.show()


### Choose same bin size as grouped_snv
def get_cumsum_throughput(df, bin_size_in_second=60):
    """
    bin_size : each bin is how many seconds
    """
    ## TODO: Max time is not always at given max. This should be fixed (if no new SNV found, copy previous row).
    df = df.sort_values('time')
    bin_size = bin_size_in_second
    # start from 1st bin or from 0 : df_both_1hr_drop_names['time'].min())-2
    # bins = range(-1, 60+1, 1)
    # labels = list(bins)
    bins = range(-bin_size, int(df['time'].max()) + bin_size, bin_size)
    labels = [f"{i}-{i + bin_size}" for i in bins[:-1]]
    df['time_bin'] = pd.cut(df['time'], bins=bins, labels=labels, right=False)

    # Group by the binned column and snv_type, then count occurrences of ALT and REF
    grouped_throughput = pd.DataFrame(df.groupby(['time_bin']).size(), columns=["read_per_bin"])
    # Calculate cumsum
    grouped_throughput['sum_reads'] = grouped_throughput['read_per_bin'].cumsum()

    return grouped_throughput, bins


def get_cumsum(df_both_1hr_drop_names, bin_size_in_second=60):
    """
    bin_size : each bin is how many seconds
    """
    # Get sample name from first entry:
    sample_name = df_both_1hr_drop_names['sample'].values[0]
    df_both_1hr_drop_names = df_both_1hr_drop_names.sort_values('time')
    bin_size = bin_size_in_second
    # Decided bins
    # bins = range(-1, 60+1, 1)
    # labels = list(bins)
    bins = range(-bin_size, int(df_both_1hr_drop_names['time'].max()) + bin_size, bin_size)
    labels = [f"{i}-{i + bin_size}" for i in bins[:-1]]
    df_both_1hr_drop_names['time_bin'] = pd.cut(df_both_1hr_drop_names['time'], bins=bins, labels=labels, right=False)

    ALT_len = df_both_1hr_drop_names[df_both_1hr_drop_names['snv_type'] == 'ALT'].shape[0]
    # In case there is no variant with ALT:
    grouped = df_both_1hr_drop_names.groupby(['time_bin', 'snv_type']).size().unstack(fill_value=0)
    if ALT_len == 0:
        print("Warning: There is no ALT in selected time period.")
        # Calculate cumsum
        grouped['ALT'] = 0
        grouped['ALT_sum'] = 0
    # Normally
    else:
        grouped['ALT_sum'] = grouped['ALT'].cumsum()
    # Group by the binned column and snv_type, then count occurrences of ALT and REF
    # Calculate cumsum

    grouped['REF_sum'] = grouped['REF'].cumsum()
    # Calculate ratio
    grouped['ratio'] = grouped['ALT_sum'] / (grouped['REF_sum'] + grouped['ALT_sum']).replace(0,
                                                                                              1)  # Replace 0 with 1 to avoid division by zero
    grouped['sum_snv'] = grouped['REF_sum'] + grouped['ALT_sum']
    grouped['sample'] = sample_name
    return grouped, bins


def plot_different_runs(df123, bins):
    df123 = df123.reset_index()
    palette_A = {
        'OVCA01': 'darkred',
        'OES01': 'darkred',
        'OES02': 'darkred',
        'GCT01': 'darkred',
        'GCT02': 'darkred',
        'HC01': 'grey',
        'HC02': 'grey',
        'HC03': 'grey',

    }

    fg = sns.catplot(x='time_bin',
                     y='ratio',
                     data=df123,
                     hue='sample',
                     aspect=4 / 2,
                     kind='point',
                     palette=palette_A, )
    locs, _ = plt.xticks()  # Get current x-tick locations
    minutes = np.array(list(bins)[1:]) / 60
    every_x = 5
    step = 0.05
    plt.yticks(ticks=np.arange(0, 1 + step, step))
    plt.xticks(locs[::every_x], labels=minutes[::every_x])  # Set x-tick locations to every 10th tick
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f1",
        "--output_pickle",
        type=str,
        help="Path to output pickle file",
        required=True,
    )
    parser.add_argument(
        "-f2",
        "--output_figure",
        type=str,
        help="Path to output figure",
        required=True,
    )

    parser.add_argument(
        "-v",
        "--vcf_name",
        type=str,
        help="vcf file sample name",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--sample",
        type=str,
        help="sample bam file used against the vcf",
        required=True,
    )
    parser.add_argument(
        "-REF",
        "--bam_ref",
        type=str,
        help="bam ref path",
        required=True,
    )
    parser.add_argument(
        "-ALT",
        "--bam_alt",
        type=str,
        help="bam alt path",
        required=True,
    )
    parser.add_argument(
        "-summary",
        "--sequencing_summary_file",
        type=str,
        help="A sequencing summary file",
        required=True,
    )
    parser.add_argument(
        "-time",
        "--time_period",
        type=int,
        help="Time period starting from sequencing",
        default=60,
    )

    args = parser.parse_args()
    time_in_minute = args.time_period
    hh = time_in_minute //60
    mm = time_in_minute % 60
    time_string = f"{hh}:{mm}:00"

    # 60 mins, each bin 1 min,
    # 180 mins, each bin 2 mins
    # 240 mins, each bin 2 mins,
    bin_size_in_second =   60
    ## Choose a period of time of which is interesting for plotting:
    sample_ref = args.vcf_name
    sample = args.sample
    # output
    output_pickle_path = args.output_pickle
    output_figure = args.output_figure
    # inputs
    sequencing_summary_file = args.sequencing_summary_file
    bam_ref = args.bam_ref
    bam_alt = args.bam_alt

    # Start
    REF_and_ALT_df = []
    sequencing_summary  = pd.read_csv(sequencing_summary_file, sep=' ', index_col=False)
    for snv_type, bam_file in zip(['REF','ALT'], [bam_ref, bam_alt]):
        snv_df = get_read_subset_with_timestamp(sequencing_summary, bam_file, sample, snv_type)
        REF_and_ALT_df.append(snv_df)
        if snv_type == "REF":
            throughput_df = get_timestamp_of_reads(sequencing_summary, sample, "all", id_list=None)
        print(f"Total {snv_type} allele SNV found (over complete sequencing summary file): ", snv_df.shape[0])
    ## Join ALT and REF:
    joined = pd.concat(REF_and_ALT_df)

    ## Method A: subsample first. Faster (2 sec)
    D0_snvs = joined[joined['timedelta'] < time_string]
    ## NEW addition: remove duplicated index to prevent duplicated index.
    D0_snvs_cumsum, D0_bins = get_cumsum(D0_snvs)
    # Get sum read count intp
    D1_reads = throughput_df[throughput_df['timedelta'] < time_string]
    D1, D1_bins= get_cumsum_throughput(D1_reads, bin_size_in_second=bin_size_in_second)
    ## Join and store the final df
    final = D1.join(D0_snvs_cumsum)
    final = final.fillna(method='ffill')
    pd.to_pickle(final, output_pickle_path)

    # # Plot alt ratio and throughput
    # # D0_snvs_cumsum may have no data on one timepoint... (but why...)
    # plot_alt_ratio(final, D1_bins, outpath)
    # # Plot cumulative count SNV
    D0_snvs = D0_snvs.reset_index(drop=True)
    plot_cum(D0_snvs, output_figure)
