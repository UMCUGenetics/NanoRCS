# NanoRCS
NanoRCS is a nanopore-based genome-wide Rolling Circle Amplification (RCA)-enhanced Consensus Sequencing method for cfDNA sequencing. 

Related publication: 
[Nanopore-based consensus sequencing enables accurate multimodal cell-free tumor DNA profiling](https://doi.org/10.1101/2024.02.16.580684).

Table of Content
==============
* [Introduction](#Introduction)
* [Repository content](#Repository-content)
  * [Process nanopore data](#Pre-processing)
  * [cfDNA Analyses](#cfDNA-Analyses)
  * [Generating figures](#Generating-figures)
* [Getting help](#Getting-help)


# Introduction
NanoRCS (pronounced "Nano rocks") is a nanopore-based genome-wide Rolling Circle Amplification (RCA)-enhanced Consensus Sequencing method for cfDNA sequencing. 

The repository was written by Li-Ting Chen and Myrthe Jager from [De Ridder lab](https://www.deridderlab.nl/) at the Center of Molecular Medicine, University Medical Center Utrecht, the Netherlands. 

You can find analysis pipelines and scripts for figures generating related to paper title: [Nanopore-based consensus sequencing enables accurate multimodal cell-free tumor DNA profiling](https://doi.org/10.1101/2024.02.16.580684).

Different parts of analyses are organized in different folders on this top folder. A README.md file could be found in each folder concerning relevant analyses. See [Repository content](#Repository content) below for content of each folder. 

The repository is still under construction for section 04-07. Please contact the authors for more details.

# Getting started:

1. Clone this repository by `git clone --recurse-submodules <address of this repo>` to include submodule ichorCNA 
and customized adaptation.
2. Install `nextflow` and `snakemake` and `R` in your environment. 
3. Further dependencies are described in each folder in this repository separately. 

# Repository content
### Pre-processing
#### 00_preprocessing_tumor_vcf
Generating a high-quality vcf file from bam files of tumor-normal pairs. These files are later used in 04_snv_tumor_informed to call SNVs in cfDNA samples.
- Related to Fig 2B-F, Fig 5A, Fig 6, Suppl Fig 2, Suppl Fig 8, Suppl Fig 10, 
#### 01_preprocessing_nanorcs
NanoRCS is composed of nanopore sequencing of RCA-enhanced cfDNA molecules.
Each nanopore read contains multiple repeats of the same cfDNA fragments. 
Consnesus algorithm is applied to generate a consensus sequence from multiple 
repeats on the same nanopore readThe settings of how we applied a consensus 
algorithm (cyclomicsseq) is detailed in this section. 

You could also find the code for raw nanorcs error rate we take the first repeat of each nanorcs read, to resemble error rate of native cfDNA sequencing. 
- Related to Fig 2A, Suppl Fig 2A, Fig 5, Suppl Fig 2, Suppl Fig 8, Suppl Fig 10, 

#### 02_preprocessing_novaseq
Illumina NovaSeq is an alternative sequencing technique we applied on majority of the cfDNA samples. A pipeline of analysing Illumina NovaSeq sequencing is provided.
- Related to Fig 2A, Fig 3C, 4D, Suppl Fig 3A-B, Suppl Fig 6
### cfDNA Analyses
#### 03_snv_error_rate
SNV error rate is derived from the bam file of 2 different sequencing techniques, with 4 different way of data processing. Results from 01,02 are required. 
- Fig 2A
#### 04_snv_tumor_informed
Tumor-informed somatic SNV detection in cfDNA samples focuses specifically on the sites where mutations have been observed in the tumor biopsy sequencing (derived from 00). 
The count of mutant alleles versus wildtype alleles were aggregated from the cfDNA molecules overlapping the mutated genomic sites, resulting in mutation fraction. 
An algorithm to derive tumor fraction in the cfDNA from the mutation fraction is implemented taken into account of variant VAF in the tumor VCF.
#### 05_cna
ichorCNA was utilized in the study to detect CNA and derive tumor fraction from CNA. Settings of ichorCNA are supplied in this folder.  
- Fig 3, 5, 6, Suppl Fig 3, 4, 8, 10
#### 06_fragmentomics
Fragmentomics is the analysis of cfDNA fragmentation patterns. We derived tumor fraction in cfDNA from fragment length distribution by NMF. 
- Fig 4, 5, 6, Suppl Fig 5,6,8,10
### 07_snv_limit_of_detection
A monte carlo simulation was used to derive the SNV limit of detection in populations of tumor patients with known tumor biopsy sequencing. 
- Suppl Fig 9

### Generating figures
#### 10_figures
All scripts for generating figures are provided per figure. A source data file can be downloaded from the supplementary information in the manuscript to generate the figures. You could also generate the source data from the related repository above.

#### output 
Default path for example output files. The files in this folder is not version tracked with git. 

# Getting help

Please contact the authors and create an issue on github to get help.


