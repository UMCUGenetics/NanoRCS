# 07_PCAWG_TF_limit_of_detection

After determining the error rate of NanoRCS, we aim to determine what is the lowest tumor fraction based on SNV observation, 
one could reliably detect at a wide-range of patients with the same cancer type. We utilize simulation to achieve this. 

This method could be applied to different sequencing methods with a known error rate. For detailed method please be referred to 
Method section "Determining lowest TF detectable with simulation on PWACG tumor patient samples", and illustration in 
Supplementary Figure 9A in the manuscript: https://doi.org/10.1101/2024.02.16.580684. 

To determine a realistic representation of a patient group of a specific cancer type, we acquired VCF files from 50 
Esophageal adenocarcinoma (EAC) patients and 100 Ovarian cancer (OVCA) patients genomic sequences from 
Pan-Cancer Analysis of Whole Genomes (PCAWG, https://dcc.icgc.org/pcawg). PCAWG contains various open-source genome sequencing of primary tumors. 

## Input  
- A known error rate per sequencing technique
- A collection of VCFs derived from patients from a specific cancer type
- Config file (see `config/` for examples): 
  - triall amount
  - technique specific sequencing throughput
  - genome length 

## Main execution steps

1. Variant calling on PCAWG patients with variant calling tool "PURPLE". Produces `${Patient_ID}.purple.somatic.postprocessed.vcf.gz`
2. Execute `scripts/vcf_filter.py` to keep single-base, PASS variants. 
3. Execute `scripts/split_patient.py` to produce `patient.pickle.gz` file to facilitate parallel execution with snakemake.
4. Execute `Snakemake-cfdetect.smk` with `config/config_YOURS.yaml` to receive detection rate per tumor fraction.
   - Output `/path/to/NanoRCS/output/processed_output/07_PCAWG_TF_limit_of_detection/combined/${TUMORTYPE}_combined.pickle.gz.csv`

## Output plotting
1. Determining reliable lowest TF detection in Fig5E and SupplFig9D: `Fig5E_SupplFig9D_v1.3.R`
2. Distribution of number of mutation in patient group `SupplFig9BC_Esophagus-cancer_PCAWG.csv`, `SupplFig9BC_Ovarian-cancer_PCAWG.csv`