###################
# This plotting script is part of the snakemake pipeline in 04_snv_tumor_informed as scripts/plot_raw_snv.py
###################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysam
import seaborn as sns
import argparse

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['font.family'] = ['Arial']

def plot_lollipop_tf(samples, tfs, tfs_err, figsize=[3,3], linewidth=4, markersize=8):
    fig, axes = plt.subplots(figsize = figsize)
    # plotting using plt.stem
    bars = axes.barh(samples,tfs,xerr=tfs_err, ) # color='white'
    axes.hlines(xmax=tfs,xmin=0, y=samples,color='grey', linewidth=4)
    axes.plot(tfs,samples, "o",color='grey', markersize = 8)
    axes.set_xlim(0, 1.0)
    axes.bar_label(bars, padding=10)
    plt.gca().invert_yaxis()
    # hide
    axes.spines[['right', 'bottom', 'top']].set_visible(False)
    # # formatting and details
    plt.xlabel('SNV-derived Tumor Fraction')


def plot_raw_snv_count(df_combined_OVCA, out_path, samples, figsize=(36, 5), plot=False, sort_by='CHROM',
                       plot_xlabel_chrom=False):
    # Old colors
    #     color_ALT = "#1456b6"
    #     color_REF = "#a9d8e0"
    #     color_ERR = "#FF2E2E"
    #     # New colors:
    color_ALT = "#A93A33"
    color_REF = "#ECECEC"
    color_ERR = "#ECECEC"

    if sort_by == 'CHROM':
        df_combined_OVCA1 = df_combined_OVCA.sort_values(['CHROM_int', 'start'], ascending=True)
        df_combined_OVCA1 = df_combined_OVCA1.reset_index(drop=True)
        ## Get chromosome position for ticks
        xtick_positions = []
        xtick_labels = []
        if plot_xlabel_chrom:
            for x in range(1, 23):
                # If no variants in one chromosome, skip.
                i = df_combined_OVCA1[df_combined_OVCA1['CHROM_int'] == x].index
                if len(i) > 0:
                    xtick_labels.append(f"chr{x}")
                    xtick_positions.append(i[0])


    elif sort_by == 'VAF':
        df_combined_OVCA1 = df_combined_OVCA.sort_values(['ALT_freq'], ascending=False)
        df_combined_OVCA1 = df_combined_OVCA1.reset_index(drop=True)

    else:
        print(
            "plot_raw_snv_count: The input dataframe is not sorted further by this function !! If want sorting, choose from CHROM or VAF.")
        #### Do not sort!!
        df_combined_OVCA1 = df_combined_OVCA.reset_index(drop=True)
    indexes = {'REF': [], 'ALT': [], 'ERR': []}
    for sample in samples:
        indexes['ALT'].append(df_combined_OVCA1[df_combined_OVCA1[sample] == 1].index.tolist())
        indexes['REF'].append(df_combined_OVCA1[df_combined_OVCA1[sample] == 0].index.tolist())
        indexes['ERR'].append(df_combined_OVCA1[df_combined_OVCA1[sample] == 2].index.tolist())

    fig, axs = plt.subplots(len(samples) + 1, 1, sharex=True, figsize=figsize)
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)
    plt.tick_params(left=False)
    # Plot each graph, and manually set the y tick values
    axs[0].spines[['right', 'bottom', 'top']].set_visible(False)
    axs[0].fill_between(range(df_combined_OVCA1.shape[0]), 0, 1, color=color_REF)
    axs[0].fill_between(range(df_combined_OVCA1.shape[0]), 0, df_combined_OVCA1.ALT_freq, color=color_ALT)
    axs[0].set_yticks(np.arange(-0.8, 1.2, 0.2))
    axs[0].set_ylim(0, 1)

    for n, sample in enumerate(samples):
        m = n + 1
        sample_index = n
        axs[m].spines[['right', 'bottom', 'top']].set_visible(False)
        axs[m].vlines(indexes['REF'][n], ymin=0.2, ymax=0.8, color=color_REF)  # linewidth=1.5
        axs[m].vlines(indexes['ALT'][n], ymin=0.2, ymax=0.8, color=color_ALT)  # linewidth=1.5
        axs[m].vlines(indexes['ERR'][n], ymin=0.2, ymax=0.8, color=color_ERR)  # linewidth=1.5
        axs[m].set_yticks(np.arange(-1, 2, 1.5), labels=["", f"{sample}"])
        axs[m].set_ylim(0, 1)

    axs[m].set_xticks(xtick_positions, labels=xtick_labels, rotation=90)
    axs[m].set_ylim(0, 1)
    plt.tight_layout()
    if sort_by == 'CHROM' and plot_xlabel_chrom:
        axs[m].xaxis.set_label_position('top')
        axs[m].set_xticks(xtick_positions, labels=xtick_labels, rotation=90)

    if plot == True:
        plt.savefig(out_path, dpi=350)
        plt.savefig(out_path + '.png', dpi=150)
    # plt.show()

    ## Second plot

    # adjust to 100%
    R1 = (df_combined_OVCA1.shape[0] - round(df_combined_OVCA1['ALT_freq'].mean() * df_combined_OVCA1.shape[0]))
    A1 = round(df_combined_OVCA1['ALT_freq'].mean() * df_combined_OVCA1.shape[0])

    stacked_portion = {"ALT": [], "REF": [], "ERR": [], "Total": []}
    for n in range(len(samples)):
        sum_up = 0
        for part in ['ALT', 'REF', 'ERR', 'Total']:

            if part != 'Total':
                stacked_portion[part].append(len(indexes[part][n]))
                part_count = len(indexes[part][n])
                sum_up += len(indexes[part][n])
            else:
                stacked_portion[part].append(sum_up)

    weight_counts = {
        "REF": [R1 / (R1 + A1) * 1000],
        "ALT": [A1 / (R1 + A1) * 1000],
        "ERR": [0 / (R1 + A1) * 1000],
    }

    for n in range(len(samples)):
        weight_counts['REF'].append(stacked_portion['REF'][n] / stacked_portion['Total'][n] * 1000)
        weight_counts['ALT'].append(stacked_portion['ALT'][n] / stacked_portion['Total'][n] * 1000)
        weight_counts['ERR'].append(stacked_portion['ERR'][n] / stacked_portion['Total'][n] * 1000)

    width = 0.1
    rows = ['Average VAF'] + samples
    fig, axs = plt.subplots(len(rows), 1, sharex=True, figsize=(6, 6))
    for item in range(len(rows)):
        bottom = 0
        p = axs[item].barh(rows[item], weight_counts['ALT'][item], width, left=bottom, color=color_ALT)
        bottom += weight_counts['ALT'][item]
        p = axs[item].barh(rows[item], weight_counts['REF'][item], width, left=bottom, color=color_REF)
        bottom += weight_counts['REF'][item]
        p = axs[item].barh(rows[item], weight_counts['ERR'][item], width, left=bottom, color=color_ERR)
        axs[item].text(1000 * 1.1,
                       0,
                       str(round(weight_counts['ALT'][item] / (weight_counts['ALT'][item] + weight_counts['REF'][item]),
                                 4)),
                       color='black',
                       )

    axs[item].set_xticks([0, 20, 100, 250, 500, 750, 900, 1000],
                         labels=["0", "2", "10", "25", "50", "75", "90", "100"], )
    plt.tight_layout()

    #
    if plot:
        plt.savefig(out_path + "bar_summary.pdf", dpi=350)
        plt.savefig(out_path + 'bar_summary.png', dpi=150)
    # plt.show()


def plot_all_vcf(all_variants,path):
    # set variant_id as index
    # all_variants.index = all_variants['variant_id']
    all_variants = all_variants[all_variants['CHROM'] != 'X']
    # all_variants = all_variants[all_variants['QUAL'] >=0]
    all_variants['var_name'] = all_variants['variant_id'].copy()
    all_variants.drop(columns=['variant_id'], inplace=True)
    all_variants = all_variants.drop_duplicates(keep='first')
    ax = sns.histplot(x = 'ALT_freq', data = all_variants, multiple="stack", color='darkred' )
    ax.set_title(f"{sample_name.split('_')[0]} All variants ALT allele frequency")
    plt.savefig(f"{path}/{sample_name.split('_')[0]}_ALT_allele_frequency-all_variants.png",dpi = 150)
    plt.savefig(f"{path}/{sample_name.split('_')[0]}_ALT_allele_frequency-all_variants.pdf")
    # plt.show()


def create_snv_table(sample_names, paths):
    data_rows = []
    dfs = []
    for sample,path in zip(sample_names,paths):
        a = pd.read_pickle(path)
        a.drop_duplicates(inplace=True)
        print("\nSample", sample)
        #     print("Unique rows (total SNV)", a.shape[0])
        b = a[a['query_map_qual'] >= 60]
        # Check GATK website for the analysis on filtering on mapping quality
        # https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
        b = b[b['CHROM'] != 'X']
        b['Sample'] = sample
        #     b = b[b['QUAL'] >=0]

        print("Ratio retained after filter for mapping quality", b.shape[0] / a.shape[0])
        #     print("TOTAL:", b.shape[0])
        dfs.append(b)

        err = b[b['snv_overlap_type'] == 'ERR'].shape[0]
        alt = b[b['snv_overlap_type'] == 'ALT'].shape[0]
        ref = b[b['snv_overlap_type'] == 'REF'].shape[0]

        data_rows.append([sample,
                          (alt / (alt + ref + err)),
                          alt,
                          (err / (alt + ref + err)),
                          err,
                          (ref / (alt + ref + err)),
                          ref,
                          (alt + ref + err),
                          ])

    with_filter_60 = pd.DataFrame(
        columns=["Sample", "ALT ratio", "ALT count", "ERR ratio", "ERR count", "REF ratio", "REF count", "TOTAL count"],
        data=data_rows)

    # X = dfs[0]
    # # Write as BED file.
    # ALT_write = X[X['snv_overlap_type'] == 'ALT'][['CHROM','start','end','REF','ALT','QUAL','FILTER','ALT_freq']]
    # ALT_write.to_csv( write_bed, sep = '\t',index=False)

    # Combine map snv type
    df_combined = pd.concat(dfs)
    df_combined = df_combined.reset_index(drop=True)
    df_combined_OVCA = df_combined.copy()

    codes = {'REF': 0, 'ALT': 1, 'ERR': 2}
    df_combined_OVCA['snv_type'] = df_combined_OVCA['snv_overlap_type'].map(codes)
    df_pivot_1 = df_combined_OVCA.pivot_table(index=['variant_id'], columns='Sample', values='snv_type', aggfunc='max')
    df_pivot_1 = df_pivot_1.fillna('-1')
    # If there are 2 reads overlapping same position in the same sample, there SNV count would be aggregated.
    df_pivot_1.reset_index(inplace=True)
    new_df_only_variant = df_combined_OVCA.drop(
        columns=['Sample', 'snv_type', 'snv_overlap_type', 'query_base', 'query_base_qual', 'query_map_qual'])
    new_df_only_variant.drop_duplicates(inplace=True)
    ab = pd.merge(right=df_pivot_1, left=new_df_only_variant, )

    df_combined_OVCA = ab.copy()
    df_combined_OVCA['CHROM_int'] = df_combined_OVCA["CHROM"].astype('float')
    df_combined_OVCA = df_combined_OVCA.sort_values(by=['CHROM_int', 'start'], ascending=True)
    df_combined_OVCA = df_combined_OVCA.reset_index(drop=True)
    for i in range(len(sample_names)):
        df_combined_OVCA[sample_names[i]] = df_combined_OVCA[sample_names[i]].astype(int)

    df_combined_OVCA.index = df_combined_OVCA['variant_id']
    df_combined_OVCA['var_name'] = df_combined_OVCA['variant_id'].copy()
    df_combined_OVCA.drop(columns=['variant_id'], inplace=True)

    return df_combined_OVCA, with_filter_60


def plot_all_observed_vaf(sample_names, paths, fig_path):
    sns.set_palette('Paired')
    # Store
    REF_name = sample_name.split('_')[0]
    for sample in sample_names:
        a_path = f"{out_vaf_observation_txt}/ref-{REF_name}-sam-{sample}-{QUAL}_adjusted-VAF.txt"
        with open(a_path, "w") as f:
            #         f.write(f"#num_measurements={df_combined[df_combined[sample]>=0].shape[0]}\n")
            df_write = df_combined_OVCA[df_combined_OVCA[sample] >= 0][['Adjusted_VAF', sample]].sort_values(
                'Adjusted_VAF')
            df_write.to_csv(a_path,
                            header=False,
                            index=False,
                            sep="\t",
                            mode='w')
        print(sample)
        print(df_combined_OVCA[df_combined_OVCA[sample] >= 0][['Adjusted_VAF', sample]].sort_values('Adjusted_VAF')[
                  sample].value_counts())

        # Plot
        subset = df_combined_OVCA[df_combined_OVCA[sample] >= 0]
        ax = sns.histplot(x='Adjusted_VAF', data=subset, hue=sample, multiple="stack", )
        # ax.legend(["ALT detected","REF detected"])
        # sns.histplot(data = x, x = 'vaf', hue = 'obs')
        ax.set_title(f"{sample}, ALT allele frequency")
        plt.savefig(
            f"{fig_path}/{sample}_{QUAL}_ALT_allele_frequency_adjusted.png",
            dpi=150)
        plt.savefig(
            f"{fig_path}/{sample}_{QUAL}_ALT_allele_frequency_adjusted.pdf")
        # plt.show()


def plot_all_vaf_in_tumor_tissue(all_variants, sample_name, path):
    ax = sns.histplot(x='ALT_freq', data=all_variants, multiple="stack",color='darkred' )
    ax.set_title(f"{sample_name.split('_')[0]} All variants ALT allele frequency")
    plt.savefig(
        f"{path}/{sample_name.split('_')[0]}_ALT_allele_frequency-all_variants.png",
        dpi=150)
    plt.savefig(
        f"{path}/{sample_name.split('_')[0]}_ALT_allele_frequency-all_variants.pdf")
    # plt.show()


def merge_df_with_all_variants(all_variants, df_combined_OVCA, sample_name):

    df_combined_OVCA_merge_all = all_variants.merge(df_combined_OVCA, how='outer')
    df_combined_OVCA_merge_all['CHROM_int'] = df_combined_OVCA_merge_all["CHROM"].astype('float')
    df_combined_OVCA_merge_all = df_combined_OVCA_merge_all.fillna('-1')
    df_combined_OVCA_merge_all.drop(columns='FILTER', inplace=True)
    df_combined_OVCA_merge_all[sample_name] = df_combined_OVCA_merge_all[sample_name].astype(int)
    assert df_combined_OVCA['var_name'].is_unique
    return df_combined_OVCA_merge_all

def plot_subset_df(df, out_path, figsize = (26, 9), method = 'min_count', count = 0):
    # Select 1200 variants to plot.
    if method == 'min_count':
        var_count = []
        for sample in sample_names:
            print(sample, df[df[sample] != -1].shape[0])
            var_count.append(df[df[sample] != -1].shape[0])
        select_count = min(var_count)
        for sample in sample_names:
            print(df[df[sample] != -1].shape[0])
            sh = df[df[sample] != -1].shape[0]
            mask = sh - select_count
            df.loc[df[df[sample] != -1].sample(n=mask, random_state=1).index, sample] = -1
            print(df[df[sample] != -1].shape)
    elif method == 'count':
        print("method count of certain number of SNV -- to be implemented/checked.")
        select_count = count
        for sample in sample_names:
            print(df[df[sample] != -1].shape[0])
            sh = df[df[sample] != -1].shape[0]
            mask = sh - select_count
            df.loc[df[df[sample] != -1].sample(n=mask, random_state=1).index, sample] = -1
            print(df[df[sample] != -1].shape)
    elif method == 'all':
        pass
    else:
        print('method is not feasible. Choose from count, min_count, and all.')
    plot_raw_snv_count(df, out_path, samples=sample_names, figsize=figsize, plot=True, sort_by='CHROM',
                       plot_xlabel_chrom=True, )


if __name__ == '__main__':
    ##################### Parameters #####################
    # Create the parser
    parser = argparse.ArgumentParser(description='Process BAM files.')
    # Add the positional argument for healthy bam files
    # nargs='+' allows for one or more arguments
    parser.add_argument('bams', nargs='+', help='List of healthy BAM files')
    # Add optional arguments
    # parser.add_argument('-n', '--sample_names', help='Name of the sample', required=False)
    parser.add_argument('-f', '--path_figure', help='Output path for the figure', required=False)
    parser.add_argument('-p', '--path', help='Output path for other figures', required=False)
    parser.add_argument('-v', '--all_variant_vcf_pickle_path', help='input all variants pickle', required=False)
    # Parse the arguments
    args = parser.parse_args()


    paths = args.bams
    sample_names = [x.split("/")[-1].split("_")[0] for x in paths]
    sample_name = sample_names[0]

    # Load pickle file for each samples included for plotting
    print("With filtering on mapping quality = ? (Depending on the config file")
    df_combined_OVCA, with_filter_60 = create_snv_table(sample_names, paths)
    # Load all molecules and merge with SNV detected molecules
    all_variants = pd.read_pickle(args.all_variant_vcf_pickle_path)
    # plot_all_vcf(all_variants)
    merge_df_with_all_variants(all_variants, df_combined_OVCA, sample_name)
    plot_all_vaf_in_tumor_tissue(all_variants, sample_name, args.path)

    # plot_all_observed_vaf((sample_names, paths, fig_path)
    plot_subset_df(df_combined_OVCA, args.path_figure, figsize = (26, 9), method = 'min_count')
