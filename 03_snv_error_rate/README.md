# cfDNA Analyses: 03_snv_error_rate

## Introduction
Analysis is related to: Method paragraph "cfDNA SNV detection error rate" in the paper.

Cell-free DNA of three HCs with known genome references were analysed. All reads overlapping SNP in any healthy controls is excluded. 

SNV error rate is derived from counting CIGAR differences of all remaining reads on chromosome 1-22 (in aggregated batch of 100,000 reads, or per read).

SNV error rate for 2 different sequencing techniques, with 4 different way of data processing are dervied. 

Output from analysis `01_preprocess_nanorcs`,`02_preprocess_novaseq` are required. 

# Steps
* Step 1: PBMC sequencing of 3 healthy controls underwent mapping and variant calling with sarek pipeline.
(`sarek_hc_variantcalling.sh`). Results from haplotypecaller is converted to bed files and merged between 3 HCs. 


* Step 2: Overlap BAM files with known VCF file.
  1. Modify config file to include relevant BAM file paths and BED file path.
  2. `snakemake --Snakefile Snakemake --cores all --configfiles configs/<config-template.yaml>` More snakemake related documentations please be referred to [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/).
  3. For RAW NanoRCS, we take all alignment (do not exclude secondary alignment). Use ` `snakemake --snakefile Snakefile_rawNanoRCS --configfiles configs/config-first-alignments.yaml --profile slurm --dry-run`
  `
* Step 3: Convert BAM files containing non-overlapping reads to pickle files that contains curated information per read. `error_rate_bam_to_pickle.py` , then generate data for plotting `error_rate_pickle_to_csv.py`
