
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
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


def mm_to_inches(mm):
    return mm * 0.0393701


def parse_tumor_fraction_from_ichorCNA(sample, input_dir):
    with open(f"{input_dir}/{sample}/{sample}.params.txt", "r") as f:
        a = f.readlines()
        for line in a:
            if line.startswith("Tumor Fraction"):
                tf = float(line.split(":")[-1].strip())
                print(sample, tf)
    return tf


def append_tf(names, input_dir):
    tfs = []
    for name in names:
        tfs.append(parse_tumor_fraction_from_ichorCNA(name, input_dir))
    return [round(tf,3) for tf in tfs]


def plot_lollipop_tf(samples, tfs, figsize=[3, 3], linewidth=4, markersize=8, output_path=None):
    """
    example: plot_lollipop_tf(['B4','B6','B9','B10','B11'], summary['SNV_derived_TF'], figsize=[3,3], linewidth=4, markersize=8)
    """
    fig, axes = plt.subplots(figsize=figsize)
    # plotting taking the idea of plt.stem, but with more controls
    bars = axes.barh(samples, tfs, color='white', lw=0)
    axes.hlines(xmax=tfs, xmin=0, y=samples, color="black", linewidth=2) # orange: #EC9D40" # pink: "#DD6FA3" # purple: "#9D58F4"
    axes.plot(tfs, samples, "o", color="black", markersize=4, )
    axes.set_xlim(0, 1.0)

    axes.bar_label(bars, padding=10)
    plt.gca().invert_yaxis()
    # hide
    axes.spines[['right', 'bottom', 'top']].set_visible(False)
    # # formatting and details
    plt.xlabel('SNV-derived Tumor Fraction')
    plt.tight_layout()
    if output_path != None:
        plt.savefig(output_path, dpi=300)


def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return r_value**2


if __name__ == '__main__':
    # Get complete cnv
    ### Define input and output
    input_dir = "/path/to/NanoRCS/output/processed_data/05_cna/ichorCNA_curated_solution/"
    output = "/path/to/NanoRCS/output/output_figures/Fig3C"
    output_tsv = "/path/to/NanoRCS/output/processed_output/05_cna"

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
    tfs = append_tf(rcs_samples, input_dir)
    temp = zip(rcs_samples, tfs)
    rcs_tf = pd.DataFrame(columns=['CYC_name', 'CYC_TF_ichorCNA'], data=temp)
    rcs_tf.to_csv(f"{output_tsv}_RCS.tsv",
                  sep='\t')
    rcs_tf.to_pickle(
        f"{output_tsv}_RCS.pickle.gz")
    tfs_ill = append_tf(ill_samples, input_dir)
    temp = zip(ill_samples, tfs_ill)
    ill_tf = pd.DataFrame(columns=['ILL_name', 'ILL_TF_ichorCNA'], data=temp)
    ill_tf.to_csv(f"{output_tsv}_ILL.tsv",
               sep='\t')
    ill_tf.to_pickle(
        f"{output_tsv}_ILL.pickle.gz")

    rcs_tf['CYC_name'] = rcs_tf['CYC_name'].str.rstrip('_CYC')
    ill_tf['ILL_name'] = ill_tf['ILL_name'].str.rstrip('_ILL')
    # plot
    lol_path = f"{output}_Fig3A_lollipop.pdf"
    lol_path_ill = f"{output}_SupplFig3A_lollipop.pdf"
    inc = mm_to_inches(6*1.13* 1.0705)
    plot_lollipop_tf(ill_tf['ILL_name'], ill_tf['ILL_TF_ichorCNA'], figsize=[inc*8, inc*8], output_path=lol_path_ill)
    plot_lollipop_tf(rcs_tf['CYC_name'], rcs_tf['CYC_TF_ichorCNA'], figsize=[inc * 3, inc * 8], output_path=lol_path)
    inc = mm_to_inches(6*1.13* 1.0705)

    # plot
    tf_df = rcs_tf.merge(ill_tf, left_on='CYC_name', right_on='ILL_name')
    def get_sample_type(sample):
        if sample.startswith("HC"):
            sample_type = 'healthy control'
        elif sample.startswith("OES"):
            sample_type = 'OES'
        elif sample.startswith("OVCA"):
            sample_type = 'OVCA'
        elif sample.startswith("GCT"):
            sample_type = 'GCT'

        return sample_type
    tf_df['Sample Type'] = tf_df['ILL_name'].apply(lambda x: get_sample_type(x))

    fig, ax = plt.subplots()
    sns.relplot(y ="CYC_TF_ichorCNA",
                x = "ILL_TF_ichorCNA",
                data= tf_df,
                hue='Sample Type',
                ax=ax,
                palette=sns.color_palette(["black", "#DD6FA3", "#9D58F4","#EC9D40",  ]),
               )
    plt.xlim(-0.03,1)
    plt.ylim(-0.03,1)

    slope, intercept, r_value, p_value, std_err = linregress(tf_df["ILL_TF_ichorCNA"], tf_df["CYC_TF_ichorCNA"])
    # slope, intercept = np.polyfit(tf_df["ILL_TF_ichorCNA"],tf_df["CYC_TF_ichorCNA"], 1)
    x = tf_df["ILL_TF_ichorCNA"]

    #add linear regression line to scatterplot
    plt.plot([0,1], [intercept,slope*1+intercept],color='grey')
    # math symbol between $$ signs. https://matplotlib.org/2.0.2/users/mathtext.html
    plt.text(0.7, 0.9, f"$R^2$={round(r_value**2,4)}")
    print(f'y={round(slope,3)}X+{round(intercept,3)}', )
    print("r-square=",round(r_value**2,4))

    plt.savefig(f"{output}_3C.pdf")
    plt.savefig(f"{output}_3C.png", dpi = 200)

