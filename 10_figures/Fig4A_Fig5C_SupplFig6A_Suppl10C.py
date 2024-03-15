import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['font.family'] = ['Arial']

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


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


def get_colors_list(samples):
    color_list = []
    color_dict = {'HC':'black','OVCA':'#DD6FA3','GCT':'#9D58F4','OES':'#EC9D40','EAC':'#EC9D40', '10percent':'black', '1percent':'black', '1':'black',  '2percent':'black','05percent':'black'}
    for sample in samples:
        color_list.append(color_dict[sample.split('0')[0]])
    return color_list


def plot_length(ps, samples, out_path, figsize=(3, 8), colors=False):

    ## Plot components:
    if colors:
        pass
    else:
        colors = get_colors_list(samples)
    i = 0
    fig, axes = plt.subplots(len(samples), figsize=figsize, sharex=True)
    # axes[0].axvline(169, color='black',linewidth=0.5)

    for sample, color in zip(samples, colors):
        axes[i].plot(ps.keys(), ps.loc[sample], label=sample, color=color)
        #     axes[i].legend(loc=(1.01,0),frameon=False)
        axes[i].spines['top'].set_visible(False)
        #     axes[i].spines['bottom'].set_visible(False)
        axes[i].spines['left'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        #     axes[i].set_ylim(0, 0.016)
        axes[i].set_yticks([])
        axes[i].set_ylabel(sample, rotation=0)
        # axes[i].s('')
        i += 1
    axes[i - 1].set_xlabel('cfDNA length')
    # plt.tight_layout()
    plt.savefig(out_path, dpi=300, transparent=True)
    plt.savefig(out_path + ".png", dpi=300, transparent=True)
    plt.show()


if __name__ == '__main__':
    ## Create figure 4A
    out_path = '/path/to/NanoRCS/output/output_figrues/Fig4A.pdf'
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
    ## plot:
    plot_length(df_CYC.iloc[:, :400-30], new_sample_names, out_path, figsize = (3,8))


    ## Create Suppl Fig. 6A
    out_path = '/path/to/NanoRCS/output/output_figrues/SupplFig6A.pdf'
    path_csv = f"/path/to/NanoRCS/output/processed_output/02_preprocessing_novaseq/results/stats_ecco/"
    suffix = ".filtered.sorted.bam.stats_len.txt"
    samples = ['HC01','HC02','HC03','OVCA01','OVCA02','OVCA03','OVCA04','OVCA05','OVCA06','OVCA07','GCT01','GCT02','OES01','OES02','OES03','OES04','OES05', ]
    df_ILL = get_combined_normalized_read_length(path_csv,
                                                 suffix,
                                                 samples,
                                                 smoothing_factor=1,
                                                 read_length_range=(31, 702))
    df_ILL.index = df_ILL.index.str.replace('OES', 'EAC')
    ## plot:
    plot_length(df_ILL.iloc[:, :400 - 30], new_sample_names, out_path, )


    ## Create Fig. 5C
    out_path = '/path/to/NanoRCS/output/output_figrues/Fig5C.pdf'
    path_csv = "/path/to/NanoRCS/output/processed_data/01_preprocess_nanorcs/consensus_filtered/stats/"
    suffix = "_len.txt"
    samples = ["OVCA01", "10percent", "2percent", "1percent", "05percent", "HC02"]
    colors = ["#DD6FA3", "#DD6FA3", "#DD6FA3", "#DD6FA3", "#DD6FA3", "#DD6FA3"]

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
                                                 smoothing_factor=1,
                                                 read_length_range=(31, 702))
    df_CYC.index = df_CYC.index.str.replace('OES', 'EAC')
    plot_length(df_CYC.iloc[:, :400 - 30], new_sample_names, out_path, figsize=(9, 4))

    ## Create Suppl Fig. 10C
    path_csv = "/Users/liting/00_projects/genome_wide_cyclomics_project/Figure4_ReadLength/data/cyclomics/"
    suffix = "_len.txt"
    samples = ["GCT02-B4",
               "GCT02-B6",
               "GCT02-B9",
               "GCT02-B10",
               "GCT02-B11",]
    colors = ["#9D58F4"] * 5 + ["black"] * 3

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
                                                 smoothing_factor=1,
                                                 read_length_range=(31, 702))
    df_CYC.index = df_CYC.index.str.replace('OES', 'EAC')

    out_path = "/path/to/NanoRCS/otuput/output_figures/SupplFig10C.pdf"
    plot_length(df_CYC.iloc[:, :800 - 30], new_sample_names, out_path, figsize=(4, 4))

