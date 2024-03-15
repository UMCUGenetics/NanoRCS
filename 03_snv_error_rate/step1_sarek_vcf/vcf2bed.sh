#!/bin/bash
GZVCF=$1
SAMPLE=$2
OUTPATH=$3
TMP=/hpc/compgen/users/lchen/tmp
THREADS=16
MAXMEM='16G'
OUTSNV=$OUTPATH/${SAMPLE}_snv.bed
OUTINS=$OUTPATH/${SAMPLE}_ins.bed
OUTDEL=$OUTPATH/${SAMPLE}_del.bed
export PATH="/hpc/compgen/users/lchen/00_utils/bin:$PATH"
# https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/vcf2bed.html
# command: bash bedops_vcf2bed.sh foo.vcf test ./
# expect output: test_snv.bed, test_ins.bed, test_del.bed
VCFnopath=$(basename -- "$GZVCF")
VCF="${OUTPATH}/${VCFnopath%.*}"

gunzip -c $GZVCF > ${VCF}
/hpc/compgen/users/lchen/00_utils/bin/convert2bed --sort-tmpdir=/hpc/compgen/users/lchen/tmp --max-mem=16G --snvs --input=vcf < $VCF > $OUTSNV
/hpc/compgen/users/lchen/00_utils/bin/convert2bed --sort-tmpdir=/hpc/compgen/users/lchen/tmp --max-mem=16G --deletions --input=vcf < $VCF > $OUTDEL
/hpc/compgen/users/lchen/00_utils/bin/convert2bed --sort-tmpdir=/hpc/compgen/users/lchen/tmp --max-mem=16G --insertions --input=vcf < $VCF > $OUTINS
cut -f1-3 $OUTSNV > ${OUTSNV%.*}_only_coor.bed
cut -f1-3 $OUTDEL > ${OUTDEL%.*}_only_coor.bed
cut -f1-3 $OUTINS > ${OUTINS%.*}_only_coor.bed
rm ${VCF}
