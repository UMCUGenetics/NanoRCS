import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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

def mm_to_inches(mm):
    return mm * 0.0393701

def get_logr_df(samples, input_dir):
    logr_copy_number = []
    ccnas = []
    for sample in samples:
        file = f"{input_dir}/{sample}/{sample}.cna.seg"
        name = file.split('/')[-1].split(".")[0]
        cna = pd.read_csv(file, sep="\t")
        cna['contig'] = cna.apply(lambda row: (row['chr'], row['start'], row['end']), axis=1)
        index = pd.MultiIndex.from_tuples(cna['contig'])
        x = cna.set_index(index)
        subset = x[[f'{name}.logR_Copy_Number']]
        subset.columns = [name]
        ccna = x[[f'{name}.Corrected_Copy_Number']]
        ccna.columns = [name]
        ccnas.append(ccna)
        logr_copy_number.append(subset)
    logr_copy_number_df = pd.concat(logr_copy_number, axis=1).T
    logr_copy_number_df.columns = logr_copy_number_df.columns.set_levels(logr_copy_number_df.columns.levels[0].astype(str), level=0)
    ccnas_df = pd.concat(ccnas, axis=1).T
    return logr_copy_number_df

def plot_pairs_supp(CYC, ILL, REF, sample):
    width = mm_to_inches(121)
    fig, axes = plt.subplots(1, 1, figsize=(width, width / 2))

    dfi = pd.DataFrame(ILL.loc[sample])
    dfc = pd.DataFrame(CYC.loc[sample])
    dft = pd.DataFrame(REF.loc[sample])

    # for visualization, all CNA > 10 is capped at 10. correlation is calculated with this capped value.
    dfc = pd.DataFrame(dfc[sample].apply(lambda x: 10 if x > 10 else x))
    dfi = pd.DataFrame(dfi[sample].apply(lambda x: 10 if x > 10 else x))
    dft = pd.DataFrame(dft[sample].apply(lambda x: 10 if x > 10 else x))

    plt.suptitle(
        f"{sample}\nPearson corr CYC-ILL: {round(dfc.corrwith(dfi).values[0], 3)}\nPearson corr CYC-TUMOR: {round(dfc.corrwith(dft).values[0], 3)}")

    # if showing overlap:
    axes.fill_between(range(dfc.T.values[0].shape[0]), [x + 4.4 for x in dfc.T.values[0]], 2 + 4.4, alpha=0.8,
                      linewidth=0.0, color="#DD6FA3", label='cfDNA NanoRCS')

    axes.fill_between(range(dfi.T.values[0].shape[0]), [x + 2.2 for x in dfi.T.values[0]], 2 + 2.2, alpha=0.8,
                      linewidth=0.0, color='lightcoral', label='cfDNA NovaSeq')

    axes.fill_between(range(dfi.T.values[0].shape[0]), [x for x in dft.T.values[0]], 2, alpha=0.5, linewidth=0.0,
                      color='grey', label='Tumor tissue')

    # axes.fill_between(range(dfi.T.values[0].shape[0]), dfc.T.values[0],  dfi.T.values[0], alpha = 1, color = 'red')

    axes.set_ylim(0, 15)
    prev = None
    # Set x labels
    xtick_pos = []
    xtick_label = []
    last_idx = 0
    for idx, key in enumerate(dfi.index):
        (contig, start, end) = key

        # Clean up contig label:
        contig = contig.replace('chr', '')
        if prev is not None and prev != contig:
            axes.axvline(idx - 0.5, c='k', lw=0.1, zorder=10)
            xtick_pos.append((idx + last_idx) / 2)
            xtick_label.append(prev)
            last_idx = idx
            # print(idx)
        prev = contig

    # Plot last tick..
    xtick_pos.append((idx + last_idx) / 2)
    xtick_label.append(contig)
    axes.set_xticks(xtick_pos, xtick_label)
    axes.legend(loc='upper right')
    # Y axis
    # axes.yaxis.tick_right()
    axes.spines['left'].set_color('grey')
    axes.tick_params(axis='y', colors='grey')
    axes.set_ylabel('CN Tumor tissue', color='grey')

    def add_shift(x):
        return x + 2.2
    def min_shift(x):
        return x - 2.2
    secax_y = axes.secondary_yaxis(
        'right', functions=(min_shift, add_shift), color='lightcoral')
    axes.spines['right'].set_color('lightcoral')

    secax_y.set_ylabel('CN cfDNA NovaSeq')
    def add_shift_2(x):
        return x + 4.4
    def min_shift_2(x):
        return x - 4.4
    # use of a float for the position:
    secax_y2 = axes.secondary_yaxis(
        1.12, functions=(min_shift_2, add_shift_2), color="#DD6FA3")
    secax_y2.set_ylabel('CN cfDNA CyclomicsSeq')
    plt.tight_layout()
    plt.savefig(f"../Figures/figure3b_{sample}_ILL_CYC_REF_supp.pdf", dpi=300)
    plt.savefig(f"../Figures/figure3b_{sample}_ILL_CYC_REF_supp.png", dpi=300)

    return 0

def plot_pairs(CYC, ILL, REF, sample):
    width = mm_to_inches(121)
    fig, axes = plt.subplots(1, 1, figsize=(width, width / 2))

    dfi = pd.DataFrame(ILL.loc[sample])
    dfc = pd.DataFrame(CYC.loc[sample])
    dft = pd.DataFrame(REF.loc[sample])

    # for visualization, all CNA > 10 is capped at 10. correlation is calculated with this capped value.
    dfc = pd.DataFrame(dfc[sample].apply(lambda x: 10 if x > 10 else x))
    dfi = pd.DataFrame(dfi[sample].apply(lambda x: 10 if x > 10 else x))
    dft = pd.DataFrame(dft[sample].apply(lambda x: 10 if x > 10 else x))

    plt.suptitle(
        f"{sample}\nPearson corr CYC-ILL: {round(dfc.corrwith(dfi).values[0], 3)}\nPearson corr CYC-TUMOR: {round(dfc.corrwith(dft).values[0], 3)}")
    #     # If showing offset:

    # if showing overlap:
    axes.fill_between(range(dfc.T.values[0].shape[0]), [x + 4.4 for x in dfc.T.values[0]], 2 + 4.4, alpha=0.8,
                      linewidth=0.0, color="#DD6FA3", label='cfDNA CyclomicsSeq')

    axes.fill_between(range(dfi.T.values[0].shape[0]), [x + 2.2 for x in dfi.T.values[0]], 2 + 2.2, alpha=0.8,
                      linewidth=0.0, color='lightcoral', label='cfDNA NovaSeq')

    axes.fill_between(range(dfi.T.values[0].shape[0]), [x for x in dft.T.values[0]], 2, alpha=0.5, linewidth=0.0,
                      color='grey', label='Tumor tissue')

    axes.set_ylim(0, 15)
    prev = None
    xtick_pos = []
    xtick_label = []
    last_idx = 0

    for idx, key in enumerate(dfi.index):
        (contig, start, end) = key

        # Clean up contig label:
        contig = contig.replace('chr', '')
        if prev is not None and prev != contig:
            axes.axvline(idx - 0.5, c='k', lw=0.1, zorder=10)
            xtick_pos.append((idx + last_idx) / 2)
            xtick_label.append(prev)
            last_idx = idx
            # print(idx)
        prev = contig
    axes.yaxis.tick_right()

    # Plot last tick..
    xtick_pos.append((idx + last_idx) / 2)
    xtick_label.append(contig)
    axes.set_xticks(xtick_pos, xtick_label)
    plt.tight_layout()
    axes.legend(loc='upper right')
    axes.spines["left"].set_visible(False)
    axes.spines["right"].set_visible(True)
    axes.spines["top"].set_visible(False)

    ## OUTPUT 1
    plt.savefig(f"../Figures/figure3b_{sample}_ILL_CYC_REF_supp.pdf", dpi=300)
    plt.savefig(f"../Figures/figure3b_{sample}_ILL_CYC_REF_supp.png", dpi=300)

    plt.show()
    return 0



if __name__ == '__main__':
    ### Define input and output
    input_dir = "/path/to/NanoRCS/output/processed_data/05_cna/ichorCNA_curated_solution/"
    output = "/path/to/NanoRCS/output/output_figures/Fig3B"

    # sample names in order:
    ill_samples = [
        "HC01_ILL",
        "HC02_ILL",
        "HC03_ILL",
        "OVCA01_ILL",
        "OVCA02_ILL",
        "OVCA03_ILL",
        "OVCA04_ILL",
        "OVCA05_ILL",
        "OVCA06_ILL",
        "OVCA07_ILL",
        "GCT01_ILL",
        "GCT02_ILL",
        "OES01_ILL",
        "OES02_ILL",
        "OES03_ILL",
        "OES04_ILL",
        "OES05_ILL",
    ]

    cyc_samples = [
        "HC01_CYC",
        "HC02_CYC",
        "HC03_CYC",
        "OVCA01_CYC",
        "OVCA02_CYC84",
        "OVCA03_CYC",
        "OVCA04_CYC",
        "OVCA05_CYC",
        "OVCA06_CYC",
        "OVCA07_CYC",
        "GCT01_CYC",
        "GCT02_CYC",
        "OES01_CYC",
        "OES02_CYC",
        "OES03_CYC",
        "OES04_CYC",
        "OES05_CYC",
    ]

    ref_samples = ['OVCA01_REF_TR2',
                   'OES01_REF',
                   'OES02_REF',
                   'GCT01_REF',
                   'GCT02_REF'
                   ]


    REF = get_logr_df(ref_samples, input_dir)
    ILL = get_logr_df(ill_samples, input_dir)
    CYC = get_logr_df(cyc_samples, input_dir)
    REF.index = REF.index.str.rstrip('_REF')
    REF.rename(index={'OVCA01_REF_TR2': 'OVCA01'}, inplace=True)
    ILL.index = ILL.index.str.rstrip('_ILL')

    CYC.index = CYC.index.str.rstrip('_CYC')
    # Comparison of raw logR value in OVCA01 for main figure
    for sample in ['OVCA01']:
        plot_pairs(CYC, ILL, REF, sample)
    # Comparison of raw logR value in all known tumor biopsy for supp figure
    for sample in ['OVCA01','OES01','OES02','GCT01','GCT02']:
        plot_pairs_supp(CYC, ILL, REF, sample)

    # Plot correlation without healthy controls.

    CYC_no_healthy = CYC.iloc[3:]
    ILL_no_healthy = ILL.iloc[3:]
    # Change this when dataset changes
    REF_no_healthy = REF.iloc[0:]
    combined1 = pd.concat([CYC_no_healthy, ILL_no_healthy])
    combined1 = combined1.T
    heat = combined1.corr(method='pearson')
    # On rows: CYC, on columns: ILL
    heat = heat.iloc[:len(CYC_no_healthy.index), len(CYC_no_healthy.index):]
    pearson_correlation = []
    for i in range(14):
        print(heat.index[i], heat.iloc[i,i], )
        pearson_correlation.append(heat.iloc[i,i])
    print("Above 0.7", sum([x>0.7 for x in pearson_correlation]))
    print("Median", np.median(pearson_correlation))
    print("Mean", np.mean(pearson_correlation))
    inc = mm_to_inches(60)
    plt.subplots(figsize=(inc, inc))
    ax = sns.heatmap(heat, cmap='Greys', vmin=0,square=True )
    ax.set_ylabel('CyclomicsSeq')
    ax.set_xlabel('NovaSeq')
    plt.tight_layout()
    plt.title('Pearson correlation of CNA per 1mb bin')
    ## OUTPUT 2

    plt.savefig(f"{output}.pdf", dpi=300)
    plt.savefig(f"{output}.png", dpi=300)

    # print("Pearson correlation:", dfi.corrwith(dfc))








