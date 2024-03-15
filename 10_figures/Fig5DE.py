import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['font.family'] = ['Arial']

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

if __name__ == '__main__':
    # INPUT
    df_sim = pd.read_csv(
        "./source_data/Fig5E.txt",
        sep='\t')

    df_all = pd.read_csv("./source_data/combined_TF_all_CYC.tsv", sep = '\t')
    # OUTPUT
    outpath = "../Figure5DE_NEW.pdf"

    # Plot dilution
    expected_values = [0.7, 0.07, 0.014, 0.007, 0.0035, ]

    df_all = df_all.set_index('sample')
    df_dilution = df_all.loc[['OVCA01','10PEROVCA01HC02','2PEROVCA01HC02','1PEROVCA01HC02','05PEROVCA01HC02', 'HC02',],:]
    df_dilution = df_dilution[['CYC Fragmentomics-derived TF','CYC CNA-derived TF','CYC SNV-derived TF']]
    df_dilution.to_csv("../dilution_series_table.tsv", sep = '\t')
    new_values = df_dilution.T.to_numpy()

    new_groups = ['Fragment length', 'CNV', 'SNV']
    dilution_levels = [ '0%', '0.5%','1%','2%', '10%', "100%",]  # Six dilution levels
    pos = np.arange(len(dilution_levels))
    bar_width = 0.25  # Adjusting the bar width for 6 groups

    # Creating 2 subplots
    fig, axes = plt.subplots(2, 1, figsize = (8, 10))
    ax = axes[0]
    # Plot subplot 1: dilution experiment. TF inferred in each modality
    for idx, group in enumerate(new_groups):
        if idx == 0:
            ax.barh(pos + idx * bar_width, new_values[idx], bar_width, label=group, left=0.001,fill=False, edgecolor ='blue', hatch='///')
        elif idx == 1:
            ax.barh(pos + idx * bar_width, new_values[idx], bar_width, label=group, left=0.001,fill=False, edgecolor ='blue')
        else:
            ax.barh(pos + idx * bar_width, new_values[idx], bar_width, label=group, left=0.001,fill=True, color='blue')


    # Setting the y-axis labels
    ax.set_yticks(pos + bar_width * len(new_groups) / 2)
    ax.set_yticklabels(dilution_levels)

    # Setting the x-axis scale and labels
    ax.set_xscale('log')
    ax.set_xlabel('Inferred tumor fraction (log)')
    ax.set_ylabel('Admixture OVCA01 ratio')
    #plt.title('Log Scale Horizontal Bar Chart with Six Groups')

    # Adjusting x-ticks to display normal numeric values
    ax.set_xticks([0.0001, 0.001, 0.01, 0.1, 1.0])
    ax.set_xticklabels(['0.0001', '0.001', '0.01', '0.1', '1.0'])

    # Adding dashed lines for the expected values at each dilution level
    for i, ev in enumerate(expected_values):
        ymin = (pos[i] - len(new_groups) *bar_width*0.1) / len(dilution_levels)
        ymax = (pos[i] + len(new_groups) * bar_width*1.1) / len(dilution_levels)
        ax.axvline(x=ev, ymin=ymin, ymax=ymax, color='grey', linestyle='--')

    ax.set_xlim(0.001, 1.0)
    ax.legend()

    ## Plot simulation SNV:
    # Define your custom palette
    custom_palette = {"Ovarian-cancer": "#DD6FA3", "Esophagus-cancer": "#EC9D40"}
    # technique order
    technique_order = ['Minion', 'NanoRCS Minion', 'Promethion', 'NanoRCS Promethion']

    ax2 = axes[1]
    sns.boxplot(y='technique',
                   x='tf',
                   hue='cancer',
                   data=df_sim,
                   palette=custom_palette,
                   linewidth=1,
                   order=technique_order,
                   fliersize=2,
                   flierprops={"marker": ".",},
                   ax = ax2,
               )
    ax2.set_xticks([0.0001, 0.001, 0.01, 0.1, 1.0])

    # Set y-axis to log scale
    ax2.legend(bbox_to_anchor=(1,1))
    ax2.set_xlim(left=0.001, right = 1.0)
    ax2.set_xscale('log')
    plt.tight_layout()
    plt.savefig(outpath)

    plt.show()

