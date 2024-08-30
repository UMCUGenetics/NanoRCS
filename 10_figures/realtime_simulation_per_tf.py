import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
# import pandas as pd
# import numpy as np

# Set font size, style
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
sns.set(rc={'font.family': 'Arial'}, font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid': False})

def get_average_read_length(read_length_table,
                            max_read_length = 1000):
    read_lengths = pd.read_csv(read_length_table, header=1)
    nanorcs_read_length = read_lengths[read_lengths['Technique'] == 'NanoRCS']
    max_read_length = max_read_length - 31
    nanorcs_average_read_length_density = nanorcs_read_length.sum(axis=0, numeric_only=True)[:max_read_length]
    return nanorcs_average_read_length_density

def get_throughput(a_folder, every_x = 5):

    ## 1. Get column "Basecalled Reads Passed" in all csv sequence throughput files in a folder.
    # Collect dfs (TODO: trim down to only valid columns if necessary)
    dfs = []
    run_names = []
    for a_file in os.listdir(a_folder):
        if a_file.endswith('.csv'):
            print(a_file)
            df1 = pd.read_csv(a_folder + a_file)
            run_name = a_file.split('/')[-1].split('.')[0]
            dfs.append(df1)
            run_names.append(run_name)
    # Collect useful columns
    for a_file in os.listdir(a_folder):
        if a_file.endswith('.csv'):
            print(a_file)
            df1 = pd.read_csv(a_folder + a_file)
            run_name = a_file.split('/')[-1].split('.')[0]
            dfs.append(df1)
            run_names.append(run_name)

    for i, df_x in enumerate(dfs):
        if i == 0:
            df_final = df_x
        else:
            df_final = df_final.join(df_x.loc[:, df1.columns.str.contains('Basecalled Reads Passed')], how='outer',
                                     rsuffix=run_names[i])
    df_final[f'Basecalled Reads Passed_{run_names[0]}'] = df_final['Basecalled Reads Passed']
    df_final.drop(columns=['Basecalled Reads Passed'], inplace=True)

    throughput_per_x_min = df_final.loc[:, df_final.columns.str.contains('Basecalled Reads Passed')]
    throughput_per_x_min = throughput_per_x_min[::every_x].copy()
    # Sequencing time = 72 hours
    throughput_per_x_min = throughput_per_x_min.loc[: 60 * 24 * 3, :]
    mean_read_throughput_per_min = throughput_per_x_min.loc[:, throughput_per_x_min.columns.str.contains('Basecalled Reads Passed')].mean(
        axis=1)[::every_x]

    return throughput_per_x_min


def get_total_cfDNA_bases(nanorcs_average_read_length_density,
                          throughput_per_x_min,
                          max_read_length = 1000,
                          fraction_reads_form_consensus = 0.5,
                          max_interested_timeframe = 360,
                          ):

    values = np.arange(0, max_read_length-31)  # from 0 to X inclusive
    fraction_of_reads_with_length_x_prob = nanorcs_average_read_length_density / np.sum(nanorcs_average_read_length_density)
    # Sampling based on known distribution in all samples (exported from Suppl Table 4)
    nested_list_for_summary_df = []
    for name, series in throughput_per_x_min[:int((max_interested_timeframe/5)+ 1)].items():
        for j, num_samples in enumerate(series):
            # print(j*5, 'mins')
            sampled_values = np.random.choice(values,
                                              size=int(num_samples* fraction_reads_form_consensus),
                                              p=fraction_of_reads_with_length_x_prob)
            total_covered_consensus_bases_per_x_mins = sampled_values.sum()
            nested_list_for_summary_df.append(
                [j * EVERY_X_MINS, num_samples, total_covered_consensus_bases_per_x_mins, name])

        df_throughput_each = pd.DataFrame(nested_list_for_summary_df,
                                     columns=['Experiment Time (minutes)', 'Basecalled Reads Passed',
                                              'Total cfDNA bases (NanoRCS consensus)', 'Run name']
                                     )
    return df_throughput_each


if __name__ == '__main__':
    ## Variables
    amount_of_tumor_specific_snv_list = [6000, 15000]
    TFs = [0.01, 0.1, 0.2, 0.4, 0.6]
    VAF = 0.5
    ## FIXED PARAMETERS for simulation error detection
    genome_size = 3_117_275_501
    nanorcs_error_rate = 0.00082 * (1/3)
    FRACTION_READS_FORM_CONSENSUS = 0.5
    # For getting sequencing data and length data
    EVERY_X_MINS = 5
    max_interested_timeframe_in_min = 360
    max_read_length = 1000
    # Plot time start, end
    time_tuple = 0, 360
    out_path = '../../output/output_figures/SupplFig3.pdf'

    ## PATHS
   a_folder = '/path/to/sequencing_summary_files'
   read_length_table = '/path/to/suppl_table_4.csv'

    ## Get read length distribution in cfDNA NanoRCS reads
    nanorcs_average_read_length_density = get_average_read_length(read_length_table,
                                                                  max_read_length = max_read_length)
    # Visualize average read length
    # plt.plot(nanorcs_average_read_length_density)
    # plt.show()


    ## Get throughput from sequencing summary files
    throughput_per_x_min = get_throughput(a_folder)

    # Throughput QC passed reads
    # sns.lineplot(data=throughput_per_x_min, palette='Purples')
    # plt.xlabel('Experiment Time (minutes)')
    # plt.ylabel('Throughput QC passed reads (n)')
    # plt.legend(bbox_to_anchor=(1, 1))


    ## Get Total cfDNA base counts by combining throughput and read length
    df_throughput_each = get_total_cfDNA_bases(nanorcs_average_read_length_density,
                          throughput_per_x_min,
                          max_read_length = max_read_length,
                          fraction_reads_form_consensus=0.5,
                          max_interested_timeframe=max_interested_timeframe_in_min,
                                               )

    ## 4. Estimate when will discovery of tumor specific SNVs is higher than background
    df_throughput = df_throughput_each.groupby(['Experiment Time (minutes)']).mean(numeric_only=True)
    df_throughput_snv1 = df_throughput.copy()
    # Create a figure with 3 subplots in a row
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel 1: Scatterplot of 'Experiment Time (minutes)' vs. 'Total cfDNA bases (NanoRCS consensus)'
    sns.scatterplot(
        x='Experiment Time (minutes)',
        y='Total cfDNA bases (NanoRCS consensus)',
        hue='Run name',
        palette='Greys',
        data=df_throughput_each,
        ax=axes[0]
    )
    axes[0].legend()
    axes[0].set_title('cfDNA Bases vs. Time')
    axes[0].set_ylabel("Total cfDNA bases (NanoRCS consensus)")
    axes[0].set_xlabel("Experiment Time (minutes)")

    num_TFs = len(TFs)
    palette1 = sns.color_palette("Purples", num_TFs)
    palette2 = sns.color_palette("Reds", num_TFs)
    palettes = [palette1, palette2]
    for i, (snv_count, palette) in enumerate(zip(amount_of_tumor_specific_snv_list, palettes)):
        for TF in TFs:
            df_throughput_snv1[f'Tumor SNV detected at {TF}'] = np.floor(
                df_throughput_snv1['Total cfDNA bases (NanoRCS consensus)'] * (
                            snv_count / genome_size) * TF * VAF)
        df_throughput_snv1[f'Tumor SNV detected at BACKGROUND'] = np.floor(
            df_throughput_snv1['Total cfDNA bases (NanoRCS consensus)'] * (
                        snv_count / genome_size) * nanorcs_error_rate)

        for color_idx, TF in enumerate(TFs):
            sns.scatterplot(
                x='Experiment Time (minutes)',
                y=f'Tumor SNV detected at {TF}',
                data=df_throughput_snv1,
                label=f'TF = {TF}',
                color = palette[color_idx],
                ax=axes[i+1]
            )
        sns.scatterplot(x='Experiment Time (minutes)', y=f'Tumor SNV detected at BACKGROUND', data=df_throughput_snv1,
                        label=f'Background errors', color = 'grey', ax=axes[i+1])

        axes[i+1].legend(bbox_to_anchor=(1, 1))
        axes[i+1].set_xlim(time_tuple)
        axes[i+1].set_title(f'Tumor SNV = {snv_count}')
        axes[i+1].set_ylabel("Tumor SNV detected (n)")
        axes[i+1].set_xlabel("Experiment Time (minutes)")

    # Panel 2: Scatterplot for Tumor SNV detected at different TFs
    # Adjust the layout to prevent overlap
    plt.tight_layout()
    # Show the final plot
    plt.savefig(out_path)
    plt.show()





    ## 6. Get the earliest SNV timepoint

    # Filter columns that start with "Tumor SNV detected at"
    tumor_snv_cols = df_throughput.loc[:, df_throughput.columns.str.contains('Tumor SNV detected at')]

    # Replace 0 with NaN
    tumor_snv_cols = tumor_snv_cols.replace(0, np.nan)

    # Find the index of the first non-zero value (now non-NaN) in each column
    first_non_zero_idx = tumor_snv_cols.idxmax()

    # Map these indices to the corresponding experiment times
    for key, value in first_non_zero_idx.items():
        print(key)
        if np.isnan(value):
            print(f'No tumor SNV discovered in maximal timespan in {key}')
        else:
            print(f'Earliest timepoint with {key}')
            first_non_zero_times = df_throughput['Experiment Time (minutes)'].iloc[int(value)]
            print(first_non_zero_times)

