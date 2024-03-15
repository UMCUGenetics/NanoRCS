# Provide output path
OUTPATH=/path/to/NanoRCS/output/processed_data/00_preprocessed_sarek_tissue_vcf/filtered_variant_calls/

#2 Filter the vcfs
python ./NanoRCS_filter_somatic_calls.py --bedtools_path /path/to/bedtools \
--temp_path /path/to/temp \
--output_path $OUTPATH \
--mutect2_prefix /path/to/NanoRCS/output/processed_data/00_preprocessed_sarek_tissue_vcf/variant-calling/SAMPLES/annotation/mutect2/ \
--strelka_prefix /path/to/NanoRCS/output/processed_data/00_preprocessed_sarek_tissue_vcf/variant-calling/SAMPLES/annotation/strelka/ \
--haplotype_caller_prefix /path/to/NanoRCS/output/processed_data/00_preprocessed_sarek_tissue_vcf/variant-calling/SAMPLES/annotation/haplotypecaller/ \
--samples_names SAMPLE01,SAMPLE02

#3 Filter the vcfs for 'PASS' only
awk '$1 ~ /^#/ || $7 == "PASS"' ${OUTPATH}/SAMPLE01_filtered_somatic.vcf > ${OUTPATH}/SAMPLE01_filtered_PASS_somatic.vcf
awk '$1 ~ /^#/ || $7 == "PASS"' ${OUTPATH}/SAMPLE02_filtered_somatic.vcf > ${OUTPATH}/SAMPLE02_filtered_PASS_somatic.vcf