#!/bin/sh
#script to run msmc on a single unphased bam file (one individual)

BAM=$1
MEAN_COVERAGE=$2
REFERENCE=$3
ROOT=`echo $1 | cut -f1 -d'.'`

samtools mpileup -q 20 -Q 20 -C 50 -u -f $REFERENCE $BAM | bcftools view -cgI - | \
/home/kdaly/bin/msmc/msmc-tools/bamCaller.py $MEAN_COVERAGE $ROOT"_mask.bed.gz" | gzip -c > $ROOT".vcf.gz"

/home/kdaly/bin/msmc/msmc-tools/generate_multihetsep.py "--mask="$ROOT"_mask.bed.gz" $BAM

