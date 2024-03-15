#!/bin/bash
## This script remove any reads that overlap with variants from the VCF file and keep the remaining BAM reads.
## BAMEXCLUDE (*_overlap_all_exclude.bam) contains bam file with reads that does not overlap with any SNV, INS, DEL.
BAM=$1
SAMPLE=$2
BEDNAME=$3
OUTPATH=$4
SNVBED=$5
INSBED=$6
DELBED=$7
TMP=/hpc/compgen/users/lchen/tmp
THREADS=16
MAXMEM='16G'
# OUTPUT
OVERLAPSNV=$OUTPATH/${SAMPLE}_${BEDNAME}_overlap_snv.bed
OVERLAPALL=$OUTPATH/${SAMPLE}_${BEDNAME}_overlap_all.bed
BAMOVERLAP=$OUTPATH/${SAMPLE}_${BEDNAME}_overlap_snv.bam
BAMEXCLUDE=$OUTPATH/${SAMPLE}_${BEDNAME}_overlap_all_exclude.bam

# Get reads in bam file overlapping with SNV file (bed format)
bedtools intersect -wa -a $BAM  -b $SNVBED > $BAMOVERLAP
# Get the coordinates of which overlaps happened.
bedtools intersect -wb -bed -a $BAM -b $SNVBED > $OVERLAPSNV
# Get reads not overlapping with SNVs (opposite from above). (-v)
bedtools intersect -wa -v -a $BAM -b $SNVBED $5 $6 > $BAMEXCLUDE 
# Get the coordinates of which overlaps happened. (Including indels)
bedtools intersect -wb -bed -a $BAM -b $SNVBED $5 $6 > $OVERLAPALL

