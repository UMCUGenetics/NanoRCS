# Preprocessing: 00_preprocessing_sarek_tissue_vcf

In this folder we perform Sarek pipeline on tumor tissue and paired normal tissue to 
acquire somatic mutation VCF files and CNA calling of "ground truth" for cell-free DNA sequencing when available.

In step1, the configuring and commands for executing sarek mapping and sarek variant calling is provided. The formatting 
of sample config file one could refer to nextflow nf-core sarek pipeline for 
specification (https://nf-co.re/sarek/3.4.0/docs/usage).

In step2, we filter and merge two VCF files called via Sarek pipeline, namely:
- mutect2 
- strelka

Variants of single-base mutation on autosomal chromosomes are kept, except the variants that are annotated as 
"germline" or "panel of normals" in mutect2. Then, PASS variants called by strelka is intersect with mutect2 set which 
forms the final set. A wrapper `step2_filter_merge_vcf/wrapper_NanoRCS_filter_somatic_calls.sh` is provided for running 
filtering on all samples.

