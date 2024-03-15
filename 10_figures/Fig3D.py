import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


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

def get_oes_related(value):
    if value in ["CCND3",
                   "VEGFA",
                   "EGFR",
                   "CDK6",
                   "GATA4",
                   "MYC",
                   "KRAS",
                   "ERBB2",
                   "GATA6",
                   "CCND1",
                   "CCNE1",
                   "MET",
                   "CDKN2A",
                   "PTEN",
                   "SMAD4",
                   "SMARCA4",
                   "APC",
                   ]:
        x = 1
    else:
        x = 0
    return x

def get_ovca_related(value):
    if value in ["KRAS",
                 "MYC",
                 "CCNE1",
                 "TP53",
                 "NF1",
                 "CDKN2A",
                 "RB1",
                 "MAP2K4",
                 "BRCA1",
                 "BRCA2",
                 ]:
        x = 1
    else:
        x = 0
    return x

def get_gct_related(value):
    # old:  "FOXL2", "DICER1","BCL11B","PRPRC","ZFHX3","CDKN2A","RB1","MAP2K4",
    if value in ["FOXL2",
                 "TP53",
                 "TERT",
                 "DICER1"
                 ]:
        x = 1
    else:
        x = 0
    return x

def get_driver_cna(value):
    if value in ["CDKN2A",
                 "PTEN",
                 "SMAD4",
                 "SMARCA4",
                 "APC",
                 "FOXL2",
                 "ZFHX3",
                 "TP53",
                 "NF1",
                 "RB1",
                 "MAP2K4",
                 "BRCA1",
                 "BRCA2"]:
        x = 'DEL'
    elif value in ["CCND3",
                   "VEGFA",
                   "EGFR",
                   "CDK6",
                   "GATA4",
                   "MYC",
                   "KRAS",
                   "ERBB2",
                   "GATA6",
                   "CCND1",
                   "CCNE1",
                   "MET",
                   "BCL11B",
                   "DICER1",
                   "PRPRC"]:
        x = 'AMP'
    else:
        x = 'unknown'
    return x


def get_gene_coordinates(input_path, output_path):
    genes = pd.read_csv(input_path, sep="\t")
    genes["gene_length"] = genes.apply(
        lambda x: int(x["Gene end (bp)"]) - int(x["Gene start (bp)"]), axis=1
    )
    genes["start"] = genes.apply(
        lambda x: int(x["Gene start (bp)"]) // 1000000 * 1000000 + 1, axis=1
    )
    genes["end"] = genes.apply(
        lambda x: ((int(x["Gene end (bp)"]) // 1000000) + 1) * 1000000, axis=1
    )
    genes["chr"] = genes["Chromosome/scaffold name"]
    if genes["Chromosome/scaffold name"].dtype == "O":
        genes = genes.drop(
            index=genes[genes["Chromosome/scaffold name"].str.contains("PATCH")].index
        )
    genes['chr'] = genes['chr'].astype(int)
    genes['driver_CNA'] = genes['Gene name'].apply(lambda x: get_driver_cna(x))
    genes['OES'] = genes['Gene name'].apply(lambda x: get_oes_related(x))
    genes['OVCA'] = genes['Gene name'].apply(lambda x: get_ovca_related(x))
    genes['GCT'] = genes['Gene name'].apply(lambda x: get_gct_related(x))
    genes = genes.sort_values(["chr","Gene start (bp)"])

    genes.reset_index(drop=True)
    genes.to_csv(output_path, sep="\t", index=False)
    print(genes)
    print("Gene list interested.")
    return genes


def plot_gene_length(genes, output_path):
    """

    :param genes: A pandas dataframe generated from get_gene_coordinates
    :return:
    """
    genes_tmp = genes.sort_values("gene_length")
    # 7.20472 inches = 183 mm
    # 3,50394 inches = 89 mm
    fig, ax = plt.subplots(figsize=[3.5, 4])
    sns.barplot(ax=ax, data=genes_tmp, x="gene_length", y="Gene name", palette="Greys")
    max_bin = round(genes_tmp["gene_length"].max() / 10000) + 2
    plt.xticks(ticks=range(0, max_bin*10000, 10000), rotation=90)
    plt.tight_layout()
    plt.savefig(
        output_path+'.pdf',
        facecolor=None,
        edgecolor=None,
        transparent=True,
        bbox_inches='tight',
        dpi=300,
    )
    plt.savefig(output_path + ".png", dpi=200)
    return 0


def _read_ichorCNA_cna_seg(genes, sample, cna, ploidy_list):

    overlapped_bins = pd.merge(
        genes,
        cna,
        how="inner",
        left_on=["chr", "start"],
        right_on=["chr", "start"],
    )
    # logR_Copy_Number is the copy number calculated from logR
    # seg.median.logR is the logR comparing to ploidy of the genome. (CNA=3 for a ploidy = 4 genome would be -0.109)
    tb = overlapped_bins[
        [
            "Gene name",
            f"{sample}.Corrected_Call",
            f"{sample}.logR_Copy_Number",
            f"{sample}.event",
            f"{sample}.Corrected_Copy_Number",
            "chr",
            "Gene start (bp)",
            "Gene end (bp)",
            "driver_CNA",
            "OES",
            "OVCA",
            "GCT"
        ]
    ].copy()
    tb.columns = [
        "gene",
        "corrected_call",
        "logR_copy_number",
        "event",
        "corrected_cnv",
        "chr",
        "start",
        "end",
        "driver_CNA",
        "OES",
        "OVCA",
        "GCT",
    ]
    tb["sample"] = sample
    if sample.startswith("HC"):
        sample_type = "Healthy"
    elif sample.startswith("OES"):
        sample_type = "OES"
    elif sample.startswith("OVCA"):
        sample_type = "OVCA"
    elif sample.startswith("GCT"):
        sample_type = "GCT"

    else:
        print(
            "Warning: sample type should be updated in def read_ichorCNA_cna_seg line 36."
        )
        exit(1)
    tb["sample_type"] = sample_type
    tb["corrected_cnv_log"] = tb['corrected_cnv'].apply(lambda x: np.log2(x/2) )

    ploidy = ploidy_list[sample]
    tb["corrected_cnv_log_ploidy_adjusted"] = tb['corrected_cnv'].apply(lambda x: np.log2(x/ploidy) )
    tb["cnv_log"] = tb['logR_copy_number'].apply(lambda x: x - 1 )
    tb["cnv_log_ploidy_adjusted"] = tb['logR_copy_number'].apply(lambda x: x - np.log2(ploidy) )


    return tb

# Get overall ploidy per sample. Use it to adjust CNV.
def get_ploidy_per_sample(input_folder, samples):
    ploidy_list = {}
    for sample in samples:
        # Per sample path, read in cna, combine with gene table.
        with open(f"{input_folder}/{sample}/{sample}.params.txt",'r') as f:
            a = f.readlines()
        for item in a:
            if item.startswith("Ploidy:"):
                ploidy_list[sample] = float(item.split("\t")[-1].strip())
    return ploidy_list

def get_all_ichorCNA_cna_seg(input_folder, samples, genes, ploidy_list):
    result_tables = []
    for sample in samples:
        # Per sample path, read in cna, combine with gene table.
        cna = pd.read_csv(f"{input_folder}/{sample}/{sample}.cna.seg", sep="\t")
        # Drop all bins on chromosome X since the gene set we are interested does not include chromosome X.
        # Otherwise, use  pd.Categorical() for sorting the chromosome from 1 to 22 and X.
        cna = cna[cna['chr'] != 'X'].copy()
        cna['chr'] = cna['chr'].astype(int)


        tb = _read_ichorCNA_cna_seg(genes, sample, cna, ploidy_list)
        result_tables.append(tb)
    combined = pd.concat(result_tables)
    combined["driver_CNA"] = pd.Categorical(
        combined["driver_CNA"], ["AMP", "DEL"]
    )

    # Remove this part
    # combined["sample_type"] = pd.Categorical(
    #     combined["sample_type"], ["OES", "OVCA", "GCT"]
    # )
    combined = combined.sort_values(['driver_CNA',
                                     'GCT',
                                     'OVCA',
                                     'OES',
                                     'gene',
                                     'sample',
                                     'sample_type'
                                     ])

    return combined


def plot_genes(
    df, output_path, color_palette=["#EC9D40", "#DD6FA3", "#9D58F4", "black"], column_names=None, colorful=False
):
    if color_palette is not None:
        color_palette = sns.color_palette(color_palette)
    p1, ax = plt.subplots(figsize=[4, 3]) # 7.2inches = 130mm

    p1 = sns.stripplot(
        ax=ax,
        data=df,
        x=column_names[0],
        y=column_names[1],
        hue=column_names[2],
        dodge=True,
        jitter=True,
        s=4,
        marker="o",
        # alpha=0.5,
        edgecolor="grey",
        linewidth=0.5,
        palette=color_palette,  # black = #050505
    )

    p1 = sns.boxplot(
        y=df[column_names[1]],
        x=df[column_names[0]],
        data=df,
        hue=column_names[2],
        palette=color_palette,
        showfliers=False,
        linewidth=0.5,
        boxprops=dict(alpha=.3),
    )
    plt.xticks(rotation=45)
    # https://github.com/mwaskom/seaborn/issues/979
    
    # For main figure, clip or crop or cat y axis above CN=10
    # p1.set_yticks([0, 2, 5, 10])
    # ax.set_ylim([0, 11])
    # For complete all datapoints
    # p1.set_yticks([0, 2, 5, 10, 15, 20, 25, 30, 35, 40])
    # ax.set_ylim([0, 41])

    plt.gca().set_xticklabels(plt.gca().get_xticklabels(), fontstyle="italic")
    # if colorful:
    #     for i, xticklabel in enumerate(plt.gca().get_xticklabels()):
    #         if xticklabel.get_text() in ["CDKN2A", "PTEN", "SMAD4","SMARCA4","APC","FOXL2",""]:
    #             xticklabel.set_color('blue')
    #         elif xticklabel.get_text() in ["CCND3", "VEGFA", "EGFR", "CDK6", "GATA4", "MYC", "KRAS", "ERBB2", "GATA6", "CCNE1", "MET"]:
    #             xticklabel.set_color('red')
    #         else:
    #             xticklabel.set_color('black')

    #  Get legend
    plt.axhline(y=0, linewidth=0.8, c='black') # label='Copy number=2'
    handles, labels = ax.get_legend_handles_labels()

    # Legend is doubled, get only the second half for circular dots , first half for box with grey lining
    l = plt.legend(
        handles[0+int(len(handles)/2):],
        labels[0+int(len(handles)/2):],
        bbox_to_anchor=(0.9 , 1),
        loc=2,
        markerscale=0.8,
        # scatterpoints=1,
        # borderaxespad=0.0,
        frameon = False,
    )
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.tight_layout()
    plt.savefig(
        output_path+'.pdf',
        facecolor=None,
        edgecolor=None,
        transparent=True,
        bbox_inches='tight',
        dpi=300,
    )
    plt.savefig(output_path + ".png", dpi=300)
    plt.show()


if __name__ == "__main__":
    input_folder = "/path/to/NanoRCS/output/processed_data/05_cna/ichorCNA_curated_solution/"
    name_prefix = "CNA-enrichment-EAC"
    ## Gene table:
    # Export from Ensembl BioMart. Included in the repo at NanoRCS/data.
    # hg19: http://grch37.ensembl.org/biomart/martview/c9a2ff2053568844ac401da986f0e895
    # hg38: https://www.ensembl.org/biomart/martview/7c1a1a47eb32480fb94270f950826fe0
    gene_table_path = "/path/to/NanoRCS/data/grch37p13_ensembl_mart_export.txt"
    # OUTPUT
    output_path = "/path/to/NanoRCS/output/output_figures/Fig3D"
    output_table = "/path/to/NanoRCS/output/processed_output/Fig3D"
    ## List of sample names to include
    # Exclude sample with TF < 0.05 which are not reliable CNA inference by ichorCNA
    samples = [
        "OES01_CYC",
        "OES02_CYC",
        "OES03_CYC",
        # "OES04_CYC",
        "OES05_CYC",
    ]
    ploidy_list = get_ploidy_per_sample(input_folder, samples)
    print(ploidy_list)
    print(f"Generating plots for {name_prefix}")
    genes = get_gene_coordinates(
        gene_table_path, f"{output_table}/gene_enrichment_test_output_df.tsv"
    )
    combined = get_all_ichorCNA_cna_seg(input_folder, samples, genes, ploidy_list)

    ## Combined subset to only OES interesting genes
    combined = combined[combined['OES'] == 1]
    # ## Crop/Clip High CNV to plot limit:
    # combined['corrected_cnv'] = combined['corrected_cnv'].clip(lower=0, upper=15)
    # Update names for plotting:
    new_names = ["Gene", "Copy number (log2 ratio)", "Sample type",]
    new_names_ploidy_adjusted = ["Gene", "Copy number ploidy adjusted (log2 ratio)", "Sample type", ]
    new_names_logR = ["Gene",  "cnv_log", "Sample type",]
    new_names_logR_ploidy_adjusted = ["Gene",  "cnv_log_ploidy_adjusted", "Sample type",]

    combined = combined.rename(
        columns={
            "corrected_cnv_log_ploidy_adjusted": new_names_ploidy_adjusted[1],
            "corrected_cnv_log": new_names[1],
            "gene": new_names[0],
            "sample_type": new_names[2],
        },
        errors="raise",
    )
    ## Plot Gene enrichment
    output_path_figures = os.path.join(output_path + f"/{name_prefix}_OES")
    # plot amplified genes
    amp_genes = combined[combined['driver_CNA'] == 'AMP']
    del_genes = combined[combined['driver_CNA'] == 'DEL']

    plot_genes(
        combined,
        output_path_figures,
        color_palette=["#EC9D40", "#DD6FA3", "#9D58F4"], #
        column_names=new_names,
        colorful=True,
    )  # ["#f77300", "#c936a0", "#8c03fc", "black"]
    print(f"Plotting logR adjusted version at output folder: {output_path}")

    # Option of plotting ploidy adjusted CNA to correct for whole genome duplication
    plot_genes(
        combined,
        output_path_figures+'_ploidy_adjusted',
        color_palette=["#EC9D40", "#DD6FA3", "#9D58F4"], #
        column_names=new_names_ploidy_adjusted,
        colorful=True,
    )  # ["#f77300", "#c936a0", "#8c03fc", "black"]
    print(f"Plotting finished without error. See plots at output folder: {output_path}")



