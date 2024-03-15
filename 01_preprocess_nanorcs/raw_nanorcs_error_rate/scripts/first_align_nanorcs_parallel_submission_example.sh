#!/bin/bash
tempbam=$1
outbam=$2
# Activate your environment include pysam
conda activate pysam_env
python first_align_nanorcs.py -i $tempbam -o $outbam
echo $outbam done
