#!/bin/sh
#Quick script to run DoC on a bam
#reference_path input_bam

out=$(ls $2 |  cut -f1 -d'_')

samtools index -@ 24 $2

java -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $1 \
--omitDepthOutputAtEachBase \
--omitIntervalStatistics \
-I $2 \
-o "DoC-autosomes_"$out \
-nt 24
-L /home/kdaly/angsd/autosomes_ARS1.list

for i in $(ls DoC* | grep -v "sample_summary"); do rm $i; done
