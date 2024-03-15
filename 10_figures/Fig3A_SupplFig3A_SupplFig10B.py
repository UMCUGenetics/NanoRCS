import numpy as np
from singlecellmultiomics.utils import is_main_chromosome, get_contig_list_from_fasta
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pysam
import seaborn as sns
import pandas as pd


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


# Define chromsome order:
def sort_chromosome_names(l):
    chrom_values = []
    for chrom in l:
        chrom_value = None
        chrom = chrom.replace('chr','').upper()
        if chrom == 'X':
            chrom_value = 99
        elif chrom == 'Y':
            chrom_value = 100
        elif chrom == 'M' or chrom=='MT':
            chrom_value = 101
        elif chrom == 'EBV':
            chrom_value = 102
        elif chrom=='MISC_ALT_CONTIGS_SCMO':
            chrom_value=999
        else:
            try:
                chrom_value = int(chrom)
            except Exception as e:
                chrom_value = 999 + sum((ord(x) for x in chrom))
        chrom_values.append(chrom_value)

    indices = sorted(range(len(chrom_values)),key=lambda x:chrom_values[x])
    return [l[idx] for idx in indices]



def parse_tumor_fraction_from_ichorCNA(sample, input_dir):
    with open(f"{input_dir}/{sample}/{sample}.params.txt", "r") as f:
        a = f.readlines()
        for line in a:
            if line.startswith("Tumor Fraction"):
                tf = float(line.split(":")[-1].strip())
                # print(tf, type(tf))
    return tf

def smoothing_cna(ccnas_df):
    # Smoothing!!!:

    df = ccnas_df.T
    df = df.fillna(2)
    pre = df.diff()
    post = df.diff(periods=-1)
    subtract = pre.subtract(post)

    df2 = df.copy()
    ## diff != 0  , subtract == 0
    for column in df.columns:
        indx1 = np.where(pre[column] != 0)
        indx2 = np.where(subtract[column] == 0)
        indx_inter = set(indx1[0]).intersection(set(indx2[0]))
        df2[column][list(indx_inter)] = df2[column][list(indx_inter)] - pre[column][list(indx_inter)]

    return df2

class GenomicPlot():
    def __init__(self, ref_path, contigs=None, ignore_contigs=None):
        """
        Initialise genomic plot

        ref_path(str or pysam.FastaFile) : Path or handle to reference

        """

        if contigs is None:
            self.contigs = sort_chromosome_names(list(
                filter(lambda x: is_main_chromosome(x) and (ignore_contigs is None or x not in ignore_contigs),
                       get_contig_list_from_fasta(ref_path))))
        else:
            self.contigs = contigs

        # Obtain the lengths:
        if type(ref_path) is str:
            with pysam.FastaFile(ref_path) as reference:
                self.lengths = {r: l for r, l in zip(reference.references, reference.lengths) if r in self.contigs}
        else:
            self.lengths = {r: l for r, l in zip(ref_path.references, ref_path.lengths) if r in self.contigs}

        self.total_bp = sum(self.lengths.values())
        # Prune contigs with no length:
        self.contigs = [contig for contig in self.contigs if contig in self.lengths]
        #         self.contigs = [contigs_conversion(contig) for contig in self.contigs]
        #         self.lengths = {contigs_conversion(r):l for (r,l) in self.lengths.items()}
        self.contigs = sort_chromosome_names(self.contigs)

    def cn_heatmap(self, df, cell_font_size=3, max_cn=4, method='ward', cmap='bwr', yticklabels=True,
                   figsize=(15, 20), xlabel='Chromosome', ylabel='Sample', vmin=0, xtickfontsize=4, **kwargs):
        """
        Create a heatmap from a copy number matrix

        df: triple indexed dataframe with as columns ('contig', start, end ), as rows cells/samples

        cell_font_size (int): font size of the cell labels

        max_cn (int) : dataframe will be clipped to this value. (Maximum copy number shown)

        method (str) : clustering metric

        cmap (str) : colormap used

        figsize(tuple) : Size of the figure

        xlabel (str) : Label for the x-axis, by default this is Contigs

        ylabel (str) : Label for the x-axis, by default this is Cells

        **kwargs : Arguments which will be passed to seaborn.clustermap

        """

        allelic_mode = len(df.columns[0]) == 4
        if allelic_mode:
            alleles = [allele for allele in df.columns.get_level_values(0).unique() if not pd.isna(allele)]
            contigs_to_plot = [contig for contig in self.contigs if contig in set(df.columns.get_level_values(1))]
            # Resample the dataframe, drop columns with no allele assigned:
            df = df.loc[:, df.columns.isin(contigs_to_plot, level=1)][alleles].sort_index(1)

            def m(k):
                allele, contig, start, end = k
                return self.contigs.index(contig), alleles.index(allele), start

            desired_order = sorted(
                list(df.loc[:, df.columns.isin(self.contigs, level=1)][alleles].sort_index(1).columns), key=m)
            df = df[desired_order]

        else:

            # Figure out what contigs are present in the dataframe:
            contigs_to_plot = [contig for contig in self.contigs if contig in set(df.columns.get_level_values(0))]
            df = df.sort_index(1)[contigs_to_plot]
        try:
            # 20230420: Changed to heatmap to make labels on the left.
            clmap = sns.heatmap(df,
                                cmap=cmap,
                                vmax=max_cn,
                                vmin=0,
                                yticklabels=yticklabels,
                                figsize=figsize,
                                **kwargs
                                )
            ax_heatmap = clmap.ax_heatmap
        except Exception as e:
            # print(e)
            # print('Falling back on heatmap without clustering')

            fig, ax_heatmap = plt.subplots(figsize=figsize)
            clmap = sns.heatmap(df,
                                cmap=cmap,
                                vmax=max_cn,
                                vmin=vmin,
                                yticklabels=True,
                                ax=ax_heatmap,
                                **kwargs
                                )

        prev = None
        xtick_pos = []
        xtick_label = []
        last_idx = 0

        allele = None
        for idx, key in enumerate(df.columns):
            if allelic_mode:
                (allele, contig, start, end) = key
            else:
                (contig, start, end) = key

            # Clean up contig label:
            contig = contig.replace('chr', '')
            if allele is not None:
                contig = f'{contig}:{allele}'
            if prev is not None and prev != contig:
                ax_heatmap.axvline(idx - 0.5, c='k', lw=0.1, zorder=10)
                xtick_pos.append((idx + last_idx) / 2)
                xtick_label.append(prev)
                last_idx = idx
                # print(idx)
            prev = contig

        # Plot last tick..
        xtick_pos.append((idx + last_idx) / 2)
        xtick_label.append(contig)

        ax_heatmap.set_xticks(xtick_pos)
        ax_heatmap.set_xticklabels(xtick_label, rotation=0, fontsize=xtickfontsize, )
        ax_heatmap.set_xlabel(xlabel, labelpad=20)
        ax_heatmap.set_ylabel(ylabel, labelpad=20)
        ax_heatmap.yaxis.set_label_position('left')


        return clmap, (ylabel, xtick_pos, xtick_label)

    def get_relative_widths(self):
        return [self.lengths[contig] / self.total_bp for contig in self.contigs]

    def reset_axis(self, contig):
        ax = self[contig]
        ax.clear()
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlabel(contig.replace('chr', ''))
        ax.set_xlim(0, self.lengths[contig])

    def get_figure(self, figsize=(20, 1)):

        widths = self.get_relative_widths()

        gs_kw = dict(width_ratios=widths)
        figure = plt.figure(figsize=figsize)
        figure.subplots_adjust(bottom=0.25, top=0.75)

        self.gridspec = gridspec.GridSpec(1, len(widths), figure=figure, wspace=0.1, width_ratios=widths)
        self.axis = {}
        prev_ax = None
        for i, contig in enumerate(self.contigs):
            # i = i + 1 # grid spec indexes from 0

            ax = plt.subplot(self.gridspec[i], sharey=prev_ax)
            self.axis[contig] = ax
            self.reset_axis(contig)
            prev_ax = ax
        sns.despine(left=True)
        figure.canvas.draw()
        return figure

    def __getitem__(self, contig):
        return self.axis[contig]


def get_df(samples, input_dir):
    result_tables = []
    logr_copy_number = []
    ccnas = []

    for sample in samples:
        file = f"{input_dir}/{sample}/{sample}.cna.seg"
        # os.listdir(input_dir):
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
        # subset = x[[f'{name}.copy.number', f'{name}.logR', f'{name}.Corrected_Copy_Number', f'{name}.logR_Copy_Number']]
    logr_copy_number_df = pd.concat(logr_copy_number, axis=1).T
    ccnas_df = pd.concat(ccnas, axis=1).T
    return ccnas_df

def append_tf(df):
    tfs = []
    for name in df.index:
        tfs.append(parse_tumor_fraction_from_ichorCNA(name, input_dir))
    return tfs



if __name__ == '__main__':
    # Define input and output
    input_dir = "/path/to/NanoRCS/output/processed_data/05_cna/ichorCNA_curated_solution/"
    output = "/path/to/NanoRCS/output/output_figures/Fig3A"


    print("Pandas version:", pd.__version__)
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

    rcs_samples = [
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

    gct_samples = ['GCT02_N',
                   'GCT02_REF',
                   'GCT02-B4_CYC',
                   'GCT02-B6_CYC',
                   'GCT02-B9_CYC',
                   'GCT02-B10_CYC',
                   'GCT02-B11_CYC']

    ref_samples = ['OVCA01_REF_TR2',
                   'OVCA01_CYC',
                   'OVCA01_ILL',
                   'GCT01_REF',
                   'GCT01_CYC',
                   'GCT01_ILL',
                   'GCT02_REF',
                   'GCT02_CYC',
                   'GCT02_ILL',
                   'OES01_REF',
                   'OES01_CYC',
                   'OES01_ILL',
                   'OES02_REF',
                   'OES02_CYC',
                   'OES02_ILL',
                   ]


    ## TO ALTER: Change this to different samples
    # rcs_samples
    # gct_samples
    # ill_samples
    ccnas_df = get_df(rcs_samples, input_dir)

    tfs = append_tf(ccnas_df)
    ccnas_df.columns = ccnas_df.columns.set_levels(ccnas_df.columns.levels[0].astype(str), level=0)

    # Smoothing!!! If a bin is different from the previous and the later bin. However the previous and the later are the same.
    # Then we apply a smoothing for this bin to become the same value as the previous bin (thus also the same as the later bin.)
    df2 = ccnas_df.T.copy()
    df2 = smoothing_cna(ccnas_df)
    df2 = df2.T

    ## Apply TF on each samples
    df2['Tf', 'tf_1', 'tf_2'] = tfs
    new_ccnas_df = df2.apply(lambda row: (row - 2) * row['Tf', 'tf_1', 'tf_2'] + 2, axis=1)
    # This line gives future warning
    new_ccnas_df.drop(columns=['Tf'], inplace=True)

    p = GenomicPlot("/path/to/references/genome/hs37d5.fa")
    x = p.cn_heatmap(new_ccnas_df, figsize=( 7.20472, 7.20472/(3) ), )
    # 7.20472 inches = 183 mm (lw = 0.2)
    # 3,50394 inches = 89 mm (lw = 0.1)
    # 20 inches (lw = 1.5)

    plt.savefig(output+".pdf", dpi = 300)
    plt.savefig(output+".png", dpi = 300)
    plt.show()
