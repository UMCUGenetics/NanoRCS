import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# This script collect multimodal tumor fraction estimation results and export tsv file.

## Sample names:
samples = ['HC01','HC02','HC03',
           'OVCA01','OVCA02','OVCA03','OVCA04','OVCA05','OVCA06','OVCA07',
           'GCT01','GCT02',
           'EAC01','EAC02','EAC03','EAC04','EAC05']

def get_sample_type(sample):
    if sample.startswith("HC"):
        sample_type = 'Healthy control'
    elif sample.startswith("OES"):
        sample_type = 'EAC'
    elif sample.startswith("EAC"):
        sample_type = 'EAC'
    elif sample.startswith("OVCA"):
        sample_type = 'OVCA'
    elif sample.startswith("GCT"):
        sample_type = 'GCT'
    elif "percent" in sample:
        sample_type = 'dilution'
    elif "OVCA01HC02" in sample:
        sample_type = 'dilution'
    return sample_type

## Mode 1: SNV

# Check if tf is between
# 0 and 0.01
# 0 and 0.1
# 0 and 1
snv_tf_path = ("/path/to/NanoRCS/output/processed_data/04_snv_tumor_informed/infer_tf/40/TF_summary")

## Healthy control against different background sample
snv_dfs = []
for healthy in ['HC01', 'HC02', 'HC03']:
    for sample in ['OVCA01', 'GCT01', 'GCT02']:
        tf_max = 0.01
        a_df = pd.read_csv(f"{snv_tf_path}/TF_max_{tf_max}_{healthy}_{sample}_tf_summary.csv", sep='\t')
        snv_dfs.append(a_df)
for healthy in ['HC01', 'HC02', 'HC03']:
    for sample in ['OES01', 'OES02']:
        tf_max = 0.1
        a_df = pd.read_csv(f"{snv_tf_path}/TF_max_{tf_max}_{healthy}_{sample}_tf_summary.csv", sep='\t')
        snv_dfs.append(a_df)

## Cancer Samples
for sample in ['OVCA01', 'GCT01', 'GCT02', 'OES01', 'OES02']:
    tf_max = 1.0
    a_df = pd.read_csv(f"{snv_tf_path}/TF_max_{tf_max}_{sample}_{sample}_tf_summary.csv", sep='\t')
    snv_dfs.append(a_df)
for sample in ['10PEROVCA01HC02']:
    tf_max = 1.0
    a_df = pd.read_csv(f"{snv_tf_path}/TF_max_{tf_max}_{sample}_{sample}_tf_summary.csv", sep='\t')
    snv_dfs.append(a_df)
for sample in ['2PEROVCA01HC02', '1PEROVCA01HC02', '05PEROVCA01HC02', ]:
    tf_max = 0.1
    a_df = pd.read_csv(f"{snv_tf_path}/TF_max_{tf_max}_{sample}_{sample}_tf_summary.csv", sep='\t')
    snv_dfs.append(a_df)

for sample in ['GCT02-B4', 'GCT02-B6', ]:
    tf_max = 1.0
    a_df = pd.read_csv(f"{snv_tf_path}/TF_max_{tf_max}_{sample}_{sample}_tf_summary.csv", sep='\t')
    snv_dfs.append(a_df)

for sample in ['GCT02-B9', 'GCT02-B10', 'GCT02-B11', ]:
    tf_max = 0.1
    a_df = pd.read_csv(f"{snv_tf_path}/TF_max_{tf_max}_{sample}_{sample}_tf_summary.csv", sep='\t')
    snv_dfs.append(a_df)

## Mode 2: Copy number alteration
p1 = "/path/to/NanoRCS/output/processed_output/05_cna/ichorCNA_derived_TF_CYC_hg19.tsv"
ichorCNA_tf = pd.read_csv(p1, sep = '\t')
p_dil = "/path/to/NanoRCS/output/processed_output/05_cna/ichorCNA_derived_TF_CYC_hg19_dilution.tsv"
ichorCNA_tf_d = pd.read_csv(p_dil, sep = '\t')
p_GCT_time = "/path/to/NanoRCS/output/processed_output/05_cna/ichorCNA_derived_TF_CYC_hg19_GCT_timeseries.tsv"
ichorCNA_tf_time = pd.read_csv(p_GCT_time, sep = '\t')
CNA_TF = pd.concat([ichorCNA_tf, ichorCNA_tf_d, ichorCNA_tf_time])
CNA_TF['name'] = CNA_TF['CYC_name'].apply(lambda x:x.split('_')[0])
CNA_TF.index = CNA_TF['name']
CNA_TF = CNA_TF[['name','CYC_TF_ichorCNA']]
CNA_TF = CNA_TF.drop_duplicates(subset = ['name'],keep='first')
CNA_TF.reset_index(drop=True, inplace= True)

## Mode 3: Fragment length
frag2 = pd.read_csv("/path/to/NanoRCS/output/processed_output/06_fragment_length_NMF_TF/NMF_length_all_samples.tsv", sep = '\t')
frag2 = frag2.drop_duplicates(subset = ['name'],keep='first')
frag2.drop(columns = ['ichorCNA_TF'], inplace = True)
frag_TF = frag2[['name','Signature2_CYC']]
frag_TF.columns = ['name', 'fragment_length_TF']


## Export tables:
# SNV only
snv_concat = pd.concat(snv_dfs)
# SNV selected columns
df2 = snv_concat[['sample','ALT_count', 'simulate_tf_max', 'SNV_derived_TF','SNV_derived_TF_range_p0025','SNV_derived_TF_range_p0975']]

SNV_TF = df2.reset_index()
SNV_TF['Sample Type'] = SNV_TF['sample'].apply(lambda x: get_sample_type(x))
SNV_TF['sample_ref'] = SNV_TF['sample']
SNV_TF['name'] = SNV_TF['sample_ref'].apply(lambda x: x.split('_')[0])
snv_tf_final = SNV_TF[['name', 'sample_ref', 'SNV_derived_TF','SNV_derived_TF_range_p0025', 'SNV_derived_TF_range_p0975']]

## Merge all
df = frag_TF.merge(CNA_TF, )
df_all = df.merge(snv_tf_final, how='outer' )


## Plot Suppl Figure 8
new_values = subset.T.to_numpy()
new_groups = ['Fragment length', 'CNV', 'SNV']
dilution_levels = subset.index # Six dilution levels
pos = np.arange(len(dilution_levels))
bar_width = 0.25  # Adjusting the bar width for 6 groups

fig, ax = plt.subplots(figsize = (10,4))

for idx, group in enumerate(new_groups):
    if idx == 0:
        ax.bar(pos + idx * bar_width, new_values[idx], bar_width/3, label=group, bottom=0.001,fill=False, edgecolor ='blue', hatch='/////')
        ax.plot(pos + idx * bar_width, new_values[idx], "o",color='blue', markersize = 6,  )
    elif idx == 1:
        ax.bar(pos + idx * bar_width, new_values[idx], bar_width/3, label=group, bottom=0.001,fill=False, edgecolor ='blue')
        ax.plot(pos + idx * bar_width, new_values[idx], "o",color='blue', markersize = 6,  )

    else:
        ax.bar(pos + idx * bar_width, new_values[idx], bar_width/3, label=group, bottom=0.001,fill=True, color='blue')
        ax.plot(pos + idx * bar_width, new_values[idx], "o",color='blue', markersize = 6,  )

# Setting the y-axis labels
ax.set_xticks(pos + bar_width * len(new_groups) / 2 - bar_width/2)
ax.set_xticklabels(dilution_levels, rotation=45)

# Setting the x-axis scale and labels
plt.ylabel('Inferred tumor fraction')
plt.xlabel('Sample Name')
plt.legend(fontsize=7, frameon=False)
## Final: Figure 6 SuppFig8.pdf
## Final: Figure 6 SuppFig8.png
plt.tight_layout()
plt.savefig("/path/to/NanoRCS/output/output_figures/SupplFig8.pdf")
plt.savefig("/path/to/NanoRCS/output/output_figures/SupplFig8.png",dpi=300)


## Take only GCT for timeseries plot
df_GCT_time_series = df_all[df_all['name'].isin(['GCT02-B4', 'GCT02-B6', 'GCT02-B9', 'GCT02-B10', 'GCT02-B11',])]

