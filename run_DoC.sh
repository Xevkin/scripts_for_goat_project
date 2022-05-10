#!/bin/sh
#Quick script to run DoC on a bam
#reference_path input_bam

out=$(ls $2 |  cut -f1 -d'_')

samtools index -@ 24 $2

/raid_md0/Software/gatk \
 DepthOfCoverage \
-R $1 \
--omit-depth-output-at-each-base \
--omit-interval-statistics \
-I $2 \
-O "DoC-autosomes_"$out \
-L /home/kdaly/raid/angsd/autosomes_ARS1.list

for i in $(ls DoC* | grep -v "sample_summary"); do rm $i; done
