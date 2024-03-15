import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict


# For plotting double axes
RATIO = 0.341
OFFSET = 20
max_plotting = 600 - 30


def nm_to_bp(nm):
    return (nm + OFFSET) / RATIO


def bp_to_nm(bp):
    return bp * RATIO - OFFSET


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


if __name__ == '__main__':
    out_path = '/path/to/NanoRCS/output/output_figrues/Fig4F'

    samples = ['HC03', 'OVCA01', 'OVCA07', 'GCT01']
    # Get Illumina df_ILL
    path_csv = f"/path/to/NanoRCS/output/processed_output/02_preprocessing_novaseq/results/stats_ecco/"
    suffix = ".filtered.sorted.bam.stats_len.txt"
    df_ILL = get_combined_normalized_read_length(path_csv,
                                                 suffix,
                                                 samples,
                                                 smoothing_factor=1,
                                                 read_length_range=(31, 702))

    # Get NanoRCS df_CYC
    path_csv = "/path/to/NanoRCS/output/processed_data/01_preprocess_nanorcs/consensus_filtered/stats/"
    suffix = "_len.txt"
    df_CYC = get_combined_normalized_read_length(path_csv,
                                             suffix,
                                             samples,
                                             smoothing_factor = 1,
                                             read_length_range=(31,701) )

    # Get AFM data & plot
    for sample in ['HC03', 'OVCA01', 'OVCA07', 'GCT01']:
        # Read AFM
        a = pd.read_csv(
            "/path/to/NanoRCS/data/Fig4_AFM_length_nm.csv")
        afm_df = a[a[sample] > 0][sample]

        # Plot
        fig, ax = plt.subplots(figsize=(4, 2.7))
        plt.title(f"{sample} length distribution")
        # Data
        secax_x = ax.secondary_xaxis(-0.2, functions=(bp_to_nm, nm_to_bp), color='saddlebrown', label='nm')
        secax_x.set_xlabel('nm')
        raw_hist = plt.hist([nm_to_bp(x) for x in afm_df], bins=200, range=(0, max_plotting), density=True, color='tan',
                            label=f'AFM (Lbp = (Lnm +{OFFSET})/{RATIO})')

        ax.set_xlabel('bp')

        plt.tight_layout()
        for df, label, color in zip([df_CYC.T, df_ILL.T], ['NanoRCS', 'NovaSeq'], ['#0070FF', '#006b3c']):
            plt.plot(df.index.values[:max_plotting], df[sample].values[:max_plotting], label=label, color=color)
        # plt.ylim(0, 0.020)
        plt.legend()
        plt.savefig(
            f"{out_path}_{sample}_r{RATIO}_o{OFFSET}_{max_plotting}.pdf")
        plt.show()