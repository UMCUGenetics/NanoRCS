output_folder="/path/to/NanoRCS/output/processed_data/03_snv_error_rate"

nextflow run /hpc/ubec/tools/pipelines/nfcore/sarek/sarek_v3.1.2/main.nf \
-profile slurm \
--custom_config_base /hpc/ubec/tools/pipelines/nfcore/sarek/sarek_v3.1.2/conf/umcu_custom \
--genome hs37d5.GRCh37 \
--step variant_calling \
--tools freebayes,strelka,haplotypecaller,mpileup,snpeff \
--input /path/to/recalibrated.csv \
--outdir ${output_folder} \
-w ${output_folder}/work \
-resume
