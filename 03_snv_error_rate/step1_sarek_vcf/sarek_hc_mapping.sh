nextflow run /path/to/pipelines/nfcore/sarek/sarek_v3.2.3.patched/main.nf -profile slurm \
--custom_config_base /path/to/nfcore/sarek/sarek_v3.4.0/conf/umcu_custom/ \
-c custom.config \
--genome hs37d5.GRCh37 \
--input /path/to/samplesheet.csv \
--step mapping \
--split_fastq 20000000 \
--trim_fasq \
--cf_window 1000 \
--outdir ${OUTPUT} \
--account ubec \
-w ${OUTPUT}/work \
-resume

if [ 0 -eq 0 ]; then
    echo 'Nextflow done.'
    echo 'NFcore workflow completed successfully.'
    rm workflow.running
    touch workflow.done
    exit 0
else
    echo 'Nextflow failed'
    rm workflow.running
    touch workflow.failed
    exit 1
fi
EOT
else
   echo 'Workflow job not submitted, please check output folder for log files.'
fi