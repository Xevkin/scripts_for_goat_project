#!/bin/bash

#script to run msmc on a single unphased bam file (one individual)
#must be indexed

BAM=$1
MEAN_COVERAGE=$2
REFERENCE=$3
CHROM_NO=$4

MASK_DIR="/kendrick/msmc/mask/"
ROOT=`echo $1 | cut -f1 -d'.'`


for i in $(eval echo "{1..$CHROM_NO}")

	do 

	echo "chr"$i
	
	samtools view -b $1 "chr"$i > $ROOT"_chr"$i".bam"

	samtools index $ROOT"_chr"$i".bam"

	~/bin/samtools-1.2/samtools mpileup -q 20 -Q 20 -C 50 -u -f $REFERENCE $ROOT"_chr"$i".bam" | bcftools call -c -V indels | \
	/home/kdaly/bin/msmc/msmc-tools/bamCaller.py $MEAN_COVERAGE  $ROOT"_chr"$i"_mask.bed.gz" | gzip -c > $ROOT"_chr"$i"_mask.vcf.gz"

	/home/kdaly/bin/msmc/msmc-tools/generate_multihetsep.py "--mask="$ROOT"_chr"$i"_mask.bed.gz" "--mask="$MASK_DIR"chr"$i".bed" $ROOT"_chr"$i"_mask.vcf.gz"  > $ROOT"_chr"$i".msmc"

done


