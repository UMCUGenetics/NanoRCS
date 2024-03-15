import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.decomposition import non_negative_factorization
import seaborn as sns

def get_sample_type(sample):
    if sample.startswith("HC"):
        sample_type = 'Healthy control'
    elif sample.startswith("OES"):
        sample_type = 'OES'
    elif sample.startswith("OVCA"):
        sample_type = 'OVCA'
    elif sample.startswith("GCT"):
        sample_type = 'GCT'
    elif "percent" in sample:
        sample_type = 'dilution'
    elif "OVCA01HC02" in sample:
        sample_type = 'dilution'
    return sample_type


def get_normalized_smoothing(cyc_dict, sample, smoothing_factor =1, read_length_range=(31, 1000)):
    d = defaultdict(int)
    for key in range(read_length_range[0], read_length_range[1]):
        # if smoothing_factor == 1, no smoothing is performed.
        written_key = key // smoothing_factor * smoothing_factor
        # If there is no value with associated length, add 0 to the key.
        if key % smoothing_factor != 0:
            try:
                d[written_key] += cyc_dict['amount'][key]
            except KeyError:
                d[written_key] += 0
        else:
            try:
                d[written_key] += cyc_dict['amount'][key]
            except KeyError:
                d[written_key] += 0
    rl = pd.DataFrame(d, index = ['amount']).T
    rl[sample] = rl['amount'].apply(lambda x: x / rl['amount'].sum())
    return rl


def get_combined_normalized_read_length(path_csv,
                                        suffix,
                                        samples,
                                        smoothing_factor = 1,
                                        read_length_range: tuple=(31, 1000)):
    all_rl = []
    for sample in samples:
        cfdna =pd.read_csv(f"{path_csv}/{sample}{suffix}",sep = "\t",comment = "#", skiprows=11, names = ['len', 'amount', 'inward','outward','other'] )
        cfdna = cfdna[['len','amount']]
        cfdna.index = cfdna['len']
        cfdna_dict = cfdna.to_dict()
        assert smoothing_factor >= 1
        rl = get_normalized_smoothing(cfdna_dict, sample,  smoothing_factor=smoothing_factor, read_length_range=read_length_range)
        all_rl.append(rl)
    df = pd.concat(all_rl, axis = 1)
    plasma_set = df[samples].T
    return plasma_set


def plot_scatter_signature_tf(df, x_column, y_column, outpath):
    sns.relplot(y=y_column,
                x=x_column,
                data=df,
                hue='Sample Type',
                ax=ax,
                palette=sns.color_palette(["black", "#DD6FA3", "#9D58F4", "#EC9D40", ]),
                )
    plt.xlim(-0.03, 1.1)
    plt.ylim(-0.03, 1.1)

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df[x_column], df[y_column])
    # slope, intercept = np.polyfit(tf_df["ILL_TF_ichorCNA"],tf_df["CYC_TF_ichorCNA"], 1)

    # add linear regression line to scatterplot
    plt.plot([0, 1], [intercept, slope * 1 + intercept], color='grey')
    # math symbol between $$ signs. https://matplotlib.org/2.0.2/users/mathtext.html
    plt.text(0.7, 0.1, f"$R^2$={round(r_value ** 2, 4)}")
    print(f'y={round(slope, 3)}X+{round(intercept, 3)}', )
    print("r-square=", round(r_value ** 2, 4))
    plt.savefig(outpath, dpi=300)
    plt.savefig(outpath + '.png', dpi=300)
    plt.show()


def plot_lollipop_tf(samples, tfs, outpath, tfs_err=None, figsize=[3,3]):
    tfs = [round(tf,3) for tf in tfs]
    fig, axes = plt.subplots(figsize = figsize)
    # plotting using plt.stem
    if tfs_err is not None:
        bars = axes.barh(samples,tfs,xerr=tfs_err, color='white')
    else:
        bars = axes.barh(samples,tfs, color='white')
    axes.hlines(xmax=tfs,xmin=0, y=samples,color='grey', linewidth=4)
    axes.plot(tfs,samples, "o",color='grey', markersize = 8)
    axes.set_xlim(0, 2.0)
    axes.bar_label(bars, padding=10)
    plt.gca().invert_yaxis()
    # hide
    axes.spines[['right', 'bottom', 'top']].set_visible(False)
    # # formatting and details
    plt.xlabel('SNV-derived Tumor Fraction')
    plt.tight_layout()
    plt.savefig(outpath)
    plt.show()

if __name__ == '__main__':
    # Plot Figure 4C
    out_path_sig = "/path/to/NanoRCS/output/output_figures/Fig4C.pdf"
    # NMF 2-signatures decomposition from WGS ILLUMINA  https://elifesciences.org/articles/71569/figures#fig2. Included in this repo
    wgs_sigs = pd.read_csv("../data/Renaud_et_al_Fig2B_WGS_Sigs_data.tsv", sep ='\t')

    until_bp_length = 280 # 700 (original signature is until 700bp. Since NanoRCS and Illumina has
    # different capture of mono- di-nucleosomal enrichment, take only read length distribution <280bp (before starting of 2nd peak).
    until_row = until_bp_length -30
    fixed_H_raw = np.array([wgs_sigs[wgs_sigs['signature_id'] == 'Signature1']['count'].values[:until_row],wgs_sigs[wgs_sigs['signature_id'] == 'Signature2']['count'].values[:until_row]])
    row_sums = fixed_H_raw.sum(axis = 1)
    ## Normalize
    fixed_H = fixed_H_raw / row_sums[:, np.newaxis]
    print(fixed_H.sum(axis = 1))
    ## plotting
    fig, ax = plt.subplots(figsize = [8,6])
    plt.plot(range(31,until_bp_length+1), fixed_H[0], color = 'grey', label='Signature 1')
    plt.plot(range(31,until_bp_length+1), fixed_H[1], color = 'darkred', label='Signature 2')
    plt.ylabel('NMF features (bp length density)')
    plt.xlabel('bp')
    plt.legend(loc=(1.0,0.8),frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path_sig)
    plt.show()



    ## Plot Suppl Fig 5
    out_path = '/path/to/NanoRCS/output/output_figrues/SupplFig5.pdf'
    # Plot Suppl Fig 5 Step1: Acquire df_CYC

    path_csv = "/path/to/NanoRCS/output/processed_data/01_preprocess_nanorcs/consensus_filtered/stats/"
    suffix = "_len.txt"
    samples = ['HC01','HC02','HC03','OVCA01','OVCA02','OVCA03','OVCA04','OVCA05','OVCA06','OVCA07','GCT01','GCT02','OES01','OES02','OES03','OES04','OES05',]
    new_sample_names = []
    for sample in samples:
        if sample.startswith('OES'):
            sample = sample.replace('OES', 'EAC')
            new_sample_names.append(sample)
        else:
            new_sample_names.append(sample)
    df_CYC = get_combined_normalized_read_length(path_csv,
                                             suffix,
                                             samples,
                                             smoothing_factor = 1,
                                             read_length_range=(31,701) )
    df_CYC.index = df_CYC.index.str.replace('OES', 'EAC')

    # Plot Suppl Fig 5 Step 2: infer NMF R matrix
    R_raw2 = df_CYC.loc[:, :until_bp_length].to_numpy()
    print(R_raw2.sum(axis=1))
    row_sums = R_raw2.sum(axis=1)
    R_CYC = R_raw2 / row_sums[:, np.newaxis]
    print(R_CYC.sum(axis=1))
    print('R shape', R_CYC.shape)

    # W: Feature matrix
    # H: Activities matrix
    W, H, n_iter = non_negative_factorization(R_CYC,
                                              n_components=2,
                                              init='custom',
                                              random_state=0,
                                              update_H=False,
                                              H=fixed_H)
    W_df = pd.DataFrame(W, columns=wgs_sigs['signature_id'].unique(), index=samples)
    W_df['Signature2'] = W_df['Signature2'].apply(lambda x: x if x < 1 else 1.0)

    ### Note: wgs_sigs['signature_id'].unique() = [Signature1, Signature2]
    # W_df is the inferred TF from NMF.

    # Plot Suppl Fig 5 Step 3: Plot length distribution against signature 1 and signature 2
    fig, axes = plt.subplots(6, 3, figsize=[3 * 5, 3 * 5])
    samples = ['HC01','HC02','HC03','OVCA01','OVCA02','OVCA03','OVCA04','OVCA05','OVCA06','OVCA07','GCT01','GCT02','OES01','OES02','OES03','OES04','OES05', ]

    for n, (item, X) in enumerate(zip(R_CYC, samples)):
        i = n // 3
        j = n % 3

        axes[i, j].plot(range(31, until_bp_length + 1),
                        wgs_sigs[wgs_sigs['signature_id'] == 'Signature1']['count'].values[:until_row],
                        color='grey',
                        linewidth=1,
                        )

        axes[i, j].plot(range(31, until_bp_length + 1),
                        wgs_sigs[wgs_sigs['signature_id'] == 'Signature2']['count'].values[:until_row],
                        color='darkred',
                        linewidth=1,
                        )
        axes[i, j].plot(range(31, until_bp_length + 1), item, linewidth=2, alpha=0.8, color='blue', label=X)

        axes[i, j].legend()
    plt.savefig(f'{out_path}')
    plt.savefig(f'{out_path}.png', dpi=300)

    ## Plot Suppl Fig 6C
    out_path_6C = '/path/to/NanoRCS/output/output_figrues/SupplFig6C.pdf'

    # Plot Suppl Fig 6C Step 1: Get df_ILL
    path_csv = f"/path/to/NanoRCS/output/processed_output/02_preprocessing_novaseq/results/stats_ecco/"
    suffix = ".filtered.sorted.bam.stats_len.txt"
    samples = ['HC01','HC02','HC03','OVCA01','OVCA02','OVCA03','OVCA04','OVCA05','OVCA06','OVCA07','GCT01','GCT02','OES01','OES02','OES03','OES04','OES05', ]
    df_ILL = get_combined_normalized_read_length(path_csv,
                                                 suffix,
                                                 samples,
                                                 smoothing_factor=1,
                                                 read_length_range=(31, 702))
    df_ILL.index = df_ILL.index.str.replace('OES', 'EAC')

    # Plot Suppl Fig 6C Step 2: Get R_ILL
    R_raw1 = df_ILL.loc[:, :until_bp_length].to_numpy()
    print(R_raw1.sum(axis=1))
    row_sums = R_raw1.sum(axis=1)
    R_ILL = R_raw1 / row_sums[:, np.newaxis]
    print(R_ILL.sum(axis=1))
    print('R shape', R_ILL.shape)

    # W: Feature matrix
    # H: Activities matrix
    W, H, n_iter = non_negative_factorization(R_ILL,
                                              n_components=2,
                                              init='custom',
                                              random_state=0,
                                              update_H=False,
                                              H=fixed_H)
    W_df_ILL = pd.DataFrame(W, columns=wgs_sigs['signature_id'].unique(), index=samples)
    W_df_ILL['Signature2'] = W_df_ILL['Signature2'].apply(lambda x: x if x < 1 else 1.0)
    df = W_df.join(W_df_ILL, lsuffix='_RCS', rsuffix='_ILL')
    df['name'] = df.index
    df['Sample Type'] = df['name'].apply(lambda x: get_sample_type(x))

    # Plot Suppl Fig 6C Step 3: Plotting ILL samples length distribution against signature 1 and signature 2
    fig, axes = plt.subplots(6, 3, figsize=[3 * 5, 3 * 5])

    for n, (item, X) in enumerate(zip(R_ILL, samples)):
        i = n // 3
        j = n % 3

        axes[i, j].plot(range(31, until_bp_length + 1),
                        wgs_sigs[wgs_sigs['signature_id'] == 'Signature1']['count'].values[:until_row],
                        color='grey',
                        linewidth=1,
                        )

        axes[i, j].plot(range(31, until_bp_length + 1),
                        wgs_sigs[wgs_sigs['signature_id'] == 'Signature2']['count'].values[:until_row],
                        color='darkred',
                        linewidth=1,
                        )
        axes[i, j].plot(range(31, until_bp_length + 1), item, linewidth=2, alpha=0.8, color='darkgreen', label=X)
        axes[i, j].legend()

    plt.savefig(f'{out_path_6C}')
    plt.savefig(f'{out_path_6C}.png', dpi=300)
    plt.show()


    ## Plot Fig 4D
    x_column = 'Signature2_ILL'
    y_column = 'Signature2_RCS'
    outpath = "/path/to/NanoRCS/output/output_figures/Fig4D.pdf"
    plot_scatter_signature_tf(df, x_column, y_column, outpath)


    ## Plot Fig 4B
    out_path_4B = "/path/to/NanoRCS/output/output_figures/Fig4B.pdf"
    plot_lollipop_tf(df['name'], df['Signature2_CYC'],
                     out_path_4B, None,
                     figsize=(3, 10))

    ## Export length inferred tf table
    df.to_csv("/path/to/NanoRCS/output/processed_output/06_fragment_length_NMF_TF/NMF_length_all_samples.tsv", sep='\t')
