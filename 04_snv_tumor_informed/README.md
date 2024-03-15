# Tumor-informed SNV analysis

VCF files from germline and somatic variant analysis in tumor biopsies with matched germline samples 
are acquired as described in `00_preprocessing_sarek_tissue_vcf`. These somatic variants specific to the tumors 
serve as the markers in tumor-informed SNV detection. 

In this folder, you find a snakemake pipeline that takes tissue VCF files, tissue tumor purity, and cfDNA bam files. 
Generating pandas dataframes and figures as output. You can find the details below.

## Method description: Monte-carlo simulation analysis
To estimate tumor fraction from cell-free DNA mixtures, we employed a Monte Carlo simulation approach. 
Cell-free DNA mixtures are composed of fragments from both healthy cells and tumor cells. 
Because of the probabilistic nature of observing a mutant (MUT) or reference (REF) allele in 
tumor-derived cfDNA. each variant position from tumor-derived cfDNA has a probability *p* of being observed 
(adjusted for tumor purity Variant Allele Frequency, VAF); for example, if the tumor purity is 60% and the 
VAF is 50%, then there is a 30% probability of observing that particular MUT allele. 

Conversely, there is a *1-p* probability of observing a REF allele. 
For healthy-cell-derived cfDNA, the probability of observing a MUT allele is 0. 
Therefore, for a given set of MUT alleles and their associated VAFs, 
we systematically vary the percentage of healthy-cell-derived cfDNA 
from 0% to 100% in 100 discrete linear steps. For each percentage, we count the number of observed MUT alleles 
by randomly sampling each allele based on their probability. We repeat this process for 10,000 trials and collect 
the observed MUT allele frequencies per tumor fraction. The most likely TF in each sample was inferred by identifying 
where the highest percentage of simulations aligned with the observed distribution of cfDNA sources. 
The confidence interval is derived from the tumor fractions that fall at the 2.5% and 97.5% of the 
simulated distribution. If the inferred TF is lower than 5%, we repeated the same process with 100 discrete log steps 
ranging from 0% to 10% to obtain a more fine-grained TF estimation. 

For each tumor tissue, we perform the Monte-carlo simulation process for targeted cfDNA sample and 
three healthy controls cfDNA samples. 

### Prerequisite:
- known tumor purity of tumor biopsy where VCF is derived. This value could be derived from purity estimation tools such as
PURPLE, CNA calling tools such as ASCAT or pathology assessment on tissue slides. 
- VCF files of tumor tissue (`00_preprocessing_sarek_tissue_vcf`)
- BAM files of cfDNA samples
- Known error rate of sequencing method (calculated from `03_snv_error_rate`)

#### For NanoRCS data processing (with multiple cfDNA from same tumor):
Execute 
```snakemake --snakefile Snakemake.smk --configfile configs/config.yaml --cores all```
- Input files please refer to example `config.yaml` file. Including
  - healthy control bam files
  - sample bam files
  - vcf for samples above
  - Repeated timepoint or dilution sample bam files
  - vcf for repeated timepoint or dilution samples
  - Tumor purity of tumor tissue

#### For NovaSeq data processing (with one cfDNA sample frome each tumor):
Execute 
```snakemake --snakefile Snakemake_ILL.smk --configfile configs/config.yaml --cores all```
- Input files please refer to example `config_ILL.yaml` file. 

### Output
Two folders contains two steps of intermediate results
- `overlap`: Contains the reads that overlap between cfDNA BAM files and the tumor tissue MUT variant positions. 
- `infer_tf`: From the overlap reads, we estimate the tumor fraction in cfDNA samples. 
- Figure 2D can be found at `overlap/{min_qual}/realtime`
- Figure 2E can be found at `overlap/{min_qual}/figures`
- Data related to Figure 2F can be found at `figure2`
The output files in `infer_tf` contains following prefix composition:
- `{cfDNA_sample_name}_{vcf_sample_name}_{TF_estimation_range_high}`

### Environment
You need to have snakemake installed in your environment (installation suggestion see snakemake website: 
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

While launching snakemake pipeline, required dependencies are included in `envs/vis.yaml`. The pipeline will automatically 
install dependencies for you.
