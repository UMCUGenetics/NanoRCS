#!/bin/bash
INPATH=$1
INFILE_NAME=$2
GENOME=$3
BACKBONE=$4

bedtools bamtofastq -i ${INPATH}/${INFILE_NAME}.bam -fq ${INFILE_NAME}.fastq
samtools view -c -F 260 ${INFILE_NAME}.fastq

cutadapt -e 0.00064 -a $BACKBONE --info-file ./cutadapt-${INFILE_NAME}.info -o ./cutadapt-trimmed-${INFILE_NAME}.fastq ${INFILE_NAME}.fastq > summary-cutadapt-${INFILE_NAME}.txt

bwa mem ${GENOME} ./cutadapt-trimmed-${INFILE_NAME}.fastq > ./cutadapt-${INFILE_NAME}.sam
samtools sort ./cutadapt-${INFILE_NAME}.sam >  ./cutadapt-${INFILE_NAME}.sorted.bam
samtools index ./cutadapt-${INFILE_NAME}.sorted.bam
rm ./cutadapt-${INFILE_NAME}.sam

